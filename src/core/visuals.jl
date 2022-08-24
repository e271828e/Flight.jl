module Visuals

using UnPack

using Printf
using CImGui
using CImGui.CSyntax
using CImGui.CSyntax.CStatic
using CImGui.GLFWBackend
using CImGui.OpenGLBackend
using CImGui.GLFWBackend.GLFW
using CImGui.OpenGLBackend.ModernGL

export CImGuiRenderer, CImGuiStyle

################################################################################
############################# CImGuiRenderer####################################

@enum CImGuiStyle begin
    classic = 0
    dark = 1
    light = 2
end

mutable struct CImGuiRenderer
    initialized::Bool
    refresh::Integer
    style::CImGuiStyle
    window::GLFW.Window
    context::Ptr{CImGui.LibCImGui.ImGuiContext}
    function CImGuiRenderer(; refresh::Integer, style::CImGuiStyle = dark)
        renderer = new()
        renderer.refresh = refresh
        renderer.style = style
        renderer.initialized = false
        return renderer
    end
end

function init!(renderer::CImGuiRenderer)

    @unpack refresh, style = renderer

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
    GLFW.SwapInterval(refresh)

    # setup Dear ImGui context
    context = CImGui.CreateContext()

    # setup Dear ImGui style
    style === classic && CImGui.StyleColorsClassic()
    style === dark && CImGui.StyleColorsDark()
    style === light && CImGui.StyleColorsLight()

    # setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true)
    ImGui_ImplOpenGL3_Init(glsl_version)

    renderer.initialized = true
    renderer.window = window
    renderer.context = context

    return nothing

end

function shutdown!(renderer::CImGuiRenderer)

    ImGui_ImplOpenGL3_Shutdown()
    ImGui_ImplGlfw_Shutdown()
    CImGui.DestroyContext(renderer.context)
    GLFW.DestroyWindow(renderer.window)
    renderer.initialized = false

    return nothing

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

    return nothing

end


function run!(renderer::CImGuiRenderer, draw::Function, draw_args...)

    renderer.initialized || init!(renderer)

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