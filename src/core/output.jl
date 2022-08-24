module Output

using Sockets
using UnPack

using Printf
using CImGui
using CImGui.CSyntax
using CImGui.CSyntax.CStatic
using CImGui.GLFWBackend
using CImGui.OpenGLBackend
using CImGui.GLFWBackend.GLFW
using CImGui.OpenGLBackend.ModernGL

export AbstractOutputInterface
export XPInterface

abstract type AbstractOutputInterface end

init!(out::AbstractOutputInterface) = throw(MethodError(init!, (out,)))
update!(out::AbstractOutputInterface, args...) = throw(MethodError(update!, (out, args...)))
# the baseline update! methods should be defined by AircraftBase

################################################################################
############################### XPInterface ####################################

Base.@kwdef struct XPInterface <: AbstractOutputInterface
    socket::UDPSocket = UDPSocket()
    host::IPv4 = IPv4("127.0.0.1")
    port::Integer = 49009
end

function set_dref(xp::XPInterface, dref_id::AbstractString, dref_value::Real)

    #ascii() ensures ASCII data, codeunits returns a CodeUnits object, which
    #behaves similarly to a byte array. this is equivalent to b"text".
    #Vector{UInt8}(dref_id) would also work
    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        UInt8(1),
        Float32(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end

function set_dref(xp::XPInterface, dref_id::AbstractString, dref_value::AbstractVector{<:Real})

    buffer = IOBuffer()
    write(buffer,
        b"DREF\0",
        dref_id |> length |> UInt8,
        dref_id |> ascii |> codeunits,
        dref_value |> length |> UInt8,
        Vector{Float32}(dref_value))

    send(xp.socket, xp.host, xp.port, buffer.data)
end

function display_text(xp::XPInterface, txt::AbstractString, x::Integer = -1, y::Integer = -1)

    buffer = IOBuffer()
    txt_ascii = ascii(txt)
    write(buffer,
        b"TEXT\0",
        x |> Int32,
        y |> Int32,
        txt_ascii |> length |> UInt8,
        txt_ascii |> codeunits)

    send(xp.socket, xp.host, xp.port, buffer.data)
end

disable_physics!(xp::XPInterface) = set_dref(xp, "sim/operation/override/override_planepath", 1)

init!(xp::XPInterface) = disable_physics!(xp)


function set_position!(xp::XPInterface; lat = -998, lon = -998, h = -998,
                        psi = -998, theta = -998, phi = -998,
                        aircraft::Integer = 0)

    #all angles must be in degrees
    buffer = IOBuffer()
    write(buffer,
        b"POSI\0", UInt8(aircraft),
        Float64(lat), Float64(lon), Float64(h),
        Float32(theta), Float32(phi), Float32(psi),
        Float32(-998)) #last one is landing gear (?!)

    send(xp.socket, xp.host, xp.port, buffer.data)

end


################################################################################
################################ CImGui ########################################

@enum CImGuiStyle begin
    classic = 0
    dark = 1
    light = 2
end

struct CImGuiRenderer{F <: Function}
    refresh_interval::Integer
    style::CImGuiStyle
    window::GLFW.Window
    context::Ptr{CImGui.LibCImGui.ImGuiContext}
end

function CImGuiRenderer(; refresh_interval::Integer, style::CImGuiStyle = dark)

    @static if Sys.isapple()
        # OpenGL 3.2 + GLSL 150
        glsl_version = 150
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 2)
        GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE) # 3.2+ only
        GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, GL_TRUE) # required on Mac
    else
        # OpenGL 3.0 + GLSL 130
        glsl_version = 130
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 0)
        # GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE) # 3.2+ only
        # GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, GL_TRUE) # 3.0+ only
    end

    # setup GLFW error callback
    error_callback(err::GLFW.GLFWError) = @error "GLFW ERROR: code $(err.code) msg: $(err.description)"
    GLFW.SetErrorCallback(error_callback)

    # create window
    window = GLFW.CreateWindow(1280, 720, "Demo")
    @assert window != C_NULL
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(refresh_interval)

    # setup Dear ImGui context
    context = CImGui.CreateContext()

    # setup Dear ImGui style
    style === classic && CImGui.StyleColorsClassic()
    style === dark && CImGui.StyleColorsDark()
    style === light && CImGui.StyleColorsLight()

    # setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true)
    ImGui_ImplOpenGL3_Init(glsl_version)

    CImGuiRenderer(refresh_interval, style, window, context)

end

function shutdown!(renderer::CImGuiRenderer)

    ImGui_ImplOpenGL3_Shutdown()
    ImGui_ImplGlfw_Shutdown()
    CImGui.DestroyContext(renderer.context)
    GLFW.DestroyWindow(renderer.window)

end


function render!(renderer::CImGuiRenderer, draw::Function, draw_args...)

    @unpack window = renderer

    try
        # start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame()
        ImGui_ImplGlfw_NewFrame()
        CImGui.NewFrame()

        #draw the frame
        draw(draw_args...)

        CImGui.Render()
        GLFW.MakeContextCurrent(window)

        display_w, display_h = GLFW.GetFramebufferSize(window)
        glViewport(0, 0, display_w, display_h)
        glClear(GL_COLOR_BUFFER_BIT)
        ImGui_ImplOpenGL3_RenderDrawData(CImGui.GetDrawData())

        GLFW.MakeContextCurrent(window)
        GLFW.SwapBuffers(window)

    catch e

        @error "Error while updating window" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
        shutdown!(renderer)

    end

end


function run!(renderer::CImGuiRenderer, draw::Function, draw_args...)

    while !GLFW.WindowShouldClose(renderer.window)
        GLFW.PollEvents()
        render!(renderer, draw, draw_args...)
    end

    shutdown!(renderer)

end


function draw_test(value::Real = 1)

    # show a simple window that we create ourselves.
    # we use a Begin/End pair to created a named window.
    @cstatic f=Cfloat(0.0) counter=Cint(0) begin
        CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
        CImGui.Text("This is some useful text.")  # display some text

        @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a slider from 0 to 1
        CImGui.Button("Button") && (counter += 1)

        CImGui.SameLine()
        CImGui.Text("counter = $counter")
        CImGui.Text("passed value = $value")
        CImGui.Text(@sprintf("Application average %.3f ms/frame (%.1f FPS)", 1000 / CImGui.GetIO().Framerate, CImGui.GetIO().Framerate))

        CImGui.End()
    end

end


# Base.@kwdef struct CImGuiDashboard{R <: CImGuiRenderer, C <: Channel} <: AbstractOutputInterface
#     renderer::R
#     channel::C
# end


#we can also define a CImGuiInput, which would have a similar interface, but
#would mimic InputManager, and would have additional arguments in draw! so that
#the input can be modified by the variables captured in the widgets


end #module