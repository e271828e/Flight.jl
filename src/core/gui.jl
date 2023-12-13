module GUI

using UnPack
using Reexport
using StaticArrays

#these are needed by any module extending GUI draw methods
@reexport using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic
@reexport using Printf
using ImGuiGLFWBackend
using ImGuiGLFWBackend.LibGLFW #defines GLFWwindow
using ImGuiGLFWBackend.LibCImGui
using ImGuiOpenGLBackend
using ImGuiOpenGLBackend.ModernGL

export CImGuiStyle, Renderer
export dynamic_button, toggle_switch, display_bar, safe_slider, safe_input
export HSV_amber, HSV_gray, HSV_green, HSV_red

################################################################################
############################### Renderer #######################################

@enum CImGuiStyle begin
    classic = 0
    dark = 1
    light = 2
end

#if sync > 0:
#T_render = T_display * sync (where typically T_display = 16.67ms).
#sync = 1 syncs the render frame rate to the display rate (vsync)
#f sync = 0)
#uncaps the render frame rate (an independent scheduling mechanism should be used)

mutable struct Renderer
    label::String
    monitor::UInt8 #which monitor to render on when multiple monitors available
    font_size::UInt8 #will be scaled by the display's content scale
    sync::UInt8 #number of display updates per frame render
    _enabled::Bool
    _initialized::Bool
    _window::Ptr{ImGuiGLFWBackend.GLFWwindow}
    _window_ctx::ImGuiGLFWBackend.Context
    _opengl_ctx::ImGuiOpenGLBackend.Context
    _cimgui_ctx::Ptr{CImGui.LibCImGui.ImGuiContext}

    function Renderer(; label = "Renderer", monitor = 2, font_size = 16, sync = 1)
        _enabled = true
        _initialized = false
        new(label, monitor, font_size, sync, _enabled, _initialized)
    end

end

Base.propertynames(::Renderer) = (:label, :monitor, :font_size, :sync)

function Base.setproperty!(renderer::Renderer, name::Symbol, value)
    if name ∈ propertynames(renderer)
        if renderer._initialized
            println("Cannot set property $name for an initialized Renderer, ",
            "call shutdown! first")
        else
            setfield!(renderer, name, value)
        end
    else
        error("Unsupported property: $name")
    end
end

enable!(renderer::Renderer) = setfield!(renderer, :_enabled, true)

function init!(renderer::Renderer)

    @unpack label, monitor, font_size, sync, _enabled = renderer

    _enabled || return

    glfwDefaultWindowHints()
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2)
    if Sys.isapple()
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE) # 3.2+ only
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE) # required on Mac
    end

    # # setup GLFW error callback
    # error_callback(err::GLFW.GLFWError) = @error "GLFW ERROR: code $(err.code) msg: $(err.description)"
    # GLFW.SetErrorCallback(error_callback)

    n_monitors = Ref{Cint}()
    monitors = glfwGetMonitors(n_monitors)
    monitor_sel = unsafe_load(monitors, min(monitor, n_monitors[])) #
    vmode = unsafe_load(glfwGetVideoMode(monitor_sel), 1)

    if n_monitors[] > 1 #create full borderless window in selected monitor
        glfwWindowHint(GLFW_FOCUSED, GLFW_TRUE)
        glfwWindowHint(GLFW_AUTO_ICONIFY, GLFW_TRUE)
        _window = glfwCreateWindow(vmode.width, vmode.height, label, monitor_sel, C_NULL)
        # _window = glfwCreateWindow(vmode.width//2, vmode.height, label, monitor_sel, C_NULL)
    else #create non-maximized window occuppying half the screen width
        glfwWindowHint(GLFW_FOCUSED, GLFW_TRUE)
        glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE)
        _window = glfwCreateWindow(vmode.width//2, vmode.height, label, C_NULL, C_NULL)
    end

    @assert _window != C_NULL
    x_pos, y_pos = Ref{Cint}(), Ref{Cint}()
    glfwGetWindowPos(_window, x_pos, y_pos)
    glfwSetWindowPos(_window, 0, y_pos[]) #no effect on borderless window

    x_scale, y_scale = Ref{Cfloat}(), Ref{Cfloat}()
    glfwGetMonitorContentScale(monitor_sel, x_scale, y_scale)
    font_scaling = max(x_scale[], y_scale[])
    scaled_font_size = round(font_scaling * font_size)

    glfwMakeContextCurrent(_window)
    glfwSwapInterval(sync)

    # create OpenGL and GLFW context
    _window_ctx = ImGuiGLFWBackend.create_context(_window)
    _opengl_ctx = ImGuiOpenGLBackend.create_context()

    # setup Dear ImGui context
    _cimgui_ctx = CImGui.CreateContext()

    # setup Dear ImGui style
    style = dark
    style === classic && CImGui.StyleColorsClassic()
    style === dark && CImGui.StyleColorsDark()
    style === light && CImGui.StyleColorsLight()

    #enable docking and multi-viewport
    #multi-viewport disabled because here it's more trouble than it's worth it
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    # io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable

    #comment when not using docking and multiviewports
    # When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(style.Colors, CImGui.ImGuiCol_WindowBg, ImVec4(col.x, col.y, col.z, 1.0f0))
    end

    fonts_dir = joinpath(@__DIR__, "gui", "fonts")
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    @assert (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Sans Linear-Regular.ttf"), scaled_font_size) != C_NULL)
    # @show (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Mono Linear-Regular.ttf"), scaled_font_size) != C_NULL)

    # setup Platform/Renderer bindings
    ImGuiGLFWBackend.init(_window_ctx)
    ImGuiOpenGLBackend.init(_opengl_ctx)

    setfield!(renderer, :_initialized, true)
    setfield!(renderer, :_window, _window)
    setfield!(renderer, :_window_ctx, _window_ctx)
    setfield!(renderer, :_opengl_ctx, _opengl_ctx)
    setfield!(renderer, :_cimgui_ctx, _cimgui_ctx)

    return nothing

end


function render(renderer::Renderer, fdraw!::Function, fdraw_args...)

    renderer._enabled || return

    @unpack _initialized, _window, _window_ctx, _opengl_ctx = renderer

    @assert _initialized "Renderer not initialized, call init! before update!"

        glfwPollEvents() #maybe on top???

        # start the Dear ImGui frame
        ImGuiOpenGLBackend.new_frame(_opengl_ctx)
        ImGuiGLFWBackend.new_frame(_window_ctx)
        CImGui.NewFrame()

        # #draw the frame and apply user inputs to arguments
        fdraw!(fdraw_args...)

        CImGui.Render()
        glfwMakeContextCurrent(_window)

        width, height = Ref{Cint}(), Ref{Cint}() #! need helper fcn
        glfwGetFramebufferSize(_window, width, height)
        display_w = width[]
        display_h = height[]

        glViewport(0, 0, display_w, display_h)
        # glClearColor(clear_color...)
        glClear(GL_COLOR_BUFFER_BIT)
        ImGuiOpenGLBackend.render(_opengl_ctx)

        #comment when not using docking and multiviewport
        if (unsafe_load(igGetIO().ConfigFlags) & ImGuiConfigFlags_ViewportsEnable) == ImGuiConfigFlags_ViewportsEnable
            ctx_backup = glfwGetCurrentContext()
            igUpdatePlatformWindows()
            GC.@preserve _opengl_ctx igRenderPlatformWindowsDefault(C_NULL, pointer_from_objref(_opengl_ctx))
            glfwMakeContextCurrent(ctx_backup)
        end

        glfwSwapBuffers(_window)

    return nothing

end


function run(renderer::Renderer, fdraw!::Function, fdraw_args...)

    renderer._enabled || return
    renderer._initialized || init!(renderer)
    try
        @assert renderer.sync > 0 "The standalone run() must not be called "*
        "an unsynced Renderer (sync = 0). Use scheduled calls to render() instead."
        #this is because a Renderer with sync=0 does not wait for monitor sync
        #when glfwSwapBuffers is called within render(), which means its frame rate
        #is effectively uncapped. this causes issues, so the calls to render must be
        #limited in frequency by some other means

        while glfwWindowShouldClose(renderer._window) == 0
            render(renderer, fdraw!, fdraw_args...)
        end
    catch e
        @error "Error while updating window" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
    finally
        shutdown!(renderer)
    end

end


function should_close(renderer::Renderer)

    renderer._enabled || return false
    renderer._initialized ? Bool(glfwWindowShouldClose(renderer._window)) : false

end


function shutdown!(renderer::Renderer)

    renderer._enabled || return
    @assert renderer._initialized "Cannot shutdown an uninitialized renderer"

    ImGuiOpenGLBackend.shutdown(renderer._opengl_ctx)
    ImGuiGLFWBackend.shutdown(renderer._window_ctx)
    CImGui.DestroyContext(renderer._cimgui_ctx)
    glfwDestroyWindow(renderer._window)
    setfield!(renderer, :_initialized, false)

    return nothing

end

function disable!(renderer::Renderer)
    !renderer._initialized ? setfield!(renderer, :_enabled, false) : println(
        "Cannot disable an already initialized renderer, call shutdown! first")
    return nothing
end

#generic non-mutating frame draw function, to be extended by users
draw(args...; kwargs...) = nothing

#generic mutating draw function, to be extended by users
draw!(args...; kwargs...) = nothing

function draw(v::AbstractVector{<:Real}, label::String, units::String = "")

    N = length(v)
    clabels = (N <= 3 ? ("X", "Y", "Z") : Tuple(1:N))

    if CImGui.TreeNode(label)
        for i in 1:N
            CImGui.Text("$(clabels[i]): ")
            CImGui.SameLine()
            CImGui.Text(@sprintf("%.7f", v[i]))
            CImGui.SameLine()
            CImGui.Text(units)
        end
        CImGui.TreePop()
    end

end

function fdraw_test(number::Real)

    CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
        output = @cstatic f=Cfloat(0.0) begin
            CImGui.Text("I got this number: $number")  # display some text
            @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a slider from 0 to 1
        end
    CImGui.End()

end

################################################################################
################################ Helpers #######################################

const HSV_gray = (0.0, 0.0, 0.3)
const HSV_amber = (0.13, 0.6, 0.6)
const HSV_green = (0.4, 0.6, 0.6)
const HSV_red = (0.0, 0.6, 0.6)

function show_help_marker(desc::String)
    CImGui.TextDisabled("(?)")
    if CImGui.IsItemHovered()
        CImGui.BeginTooltip()
        CImGui.PushTextWrapPos(CImGui.GetFontSize() * 35.0)
        CImGui.TextUnformatted(desc)
        CImGui.PopTextWrapPos()
        CImGui.EndTooltip()
    end
end

#changes shade when hovered and pushed
function toggle_switch(label::String, hue::AbstractFloat, enabled::Bool)
    if enabled
        CImGui.PushStyleColor(CImGui.ImGuiCol_Button, CImGui.HSV(hue, 0.8, 0.8))
    else
        CImGui.PushStyleColor(CImGui.ImGuiCol_Button, CImGui.HSV(hue, 0.2, 0.2))
    end
    CImGui.Button(label)
    CImGui.PopStyleColor(1)
    enable = CImGui.IsItemActive()
    return enable
end

function dynamic_button(label::String,
                        idle_HSV::NTuple{3,Real},
                        hover_HSV::NTuple{3,Real},
                        push_HSV::NTuple{3,Real})

    idle_HSV_bnd = min.(max.(SVector{3, Float64}(idle_HSV), 0.0), 1.0)
    hover_HSV_bnd = min.(max.(SVector{3, Float64}(hover_HSV), 0.0), 1.0)
    push_HSV_bnd = min.(max.(SVector{3, Float64}(push_HSV), 0.0), 1.0)
    CImGui.PushStyleColor(CImGui.ImGuiCol_Button, CImGui.HSV(idle_HSV_bnd...))
    CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonHovered, CImGui.HSV(hover_HSV_bnd...))
    CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonActive, CImGui.HSV(push_HSV_bnd...))
    CImGui.Button(label)
    CImGui.PopStyleColor(3)
    return nothing
end

function dynamic_button(label::String,
                        idle_HSV::NTuple{3,Real},
                        Δ_hover::Real = 0.1,
                        Δ_push::Real = 0.2)

    hover_HSV = (idle_HSV[1], idle_HSV[2] + Δ_hover, idle_HSV[3] + Δ_hover)
    push_HSV = (idle_HSV[1], idle_HSV[2] + Δ_push, idle_HSV[3] + Δ_push)
    dynamic_button(label, idle_HSV, hover_HSV, push_HSV)

end

#changes shade when hovered and pushed
function dynamic_button(label::String, hue::AbstractFloat = 0.4)
    dynamic_button(label, (hue, 0.6, 0.6))
end


function display_bar(label::String, source::Real, lower_bound::Real, upper_bound::Real, size_arg = (0, 0))
    CImGui.Text(label)
    CImGui.SameLine()
    CImGui.ProgressBar((source - lower_bound)/(upper_bound - lower_bound), size_arg, "$source")
end

#the string after ## is not shown, but is part of the widget's ID, which must be
#unique to avoid conflicts. an alternative solution is to use PushID and PopID.
#See DearImGui FAQ
function safe_slider(label::String, source::AbstractFloat, lower_bound::Real, upper_bound::Real, display_format::String, show_label = false)
    ref = Ref(Cfloat(source))
    slider_label = show_label ? label : "##"*label
    CImGui.SliderFloat(slider_label, ref, lower_bound, upper_bound, display_format)
    # CImGui.SameLine()
    # show_help_marker("Ctrl+Click for keyboard input")
    return ref[]
end

function safe_input(label::String, source::AbstractFloat, step::Real, fast_step::Real, display_format::String, show_label = false)
    ref = Ref(Cdouble(source))
    input_label = show_label ? label : "##"*label
    CImGui.InputDouble(input_label, ref, step, fast_step, display_format)
    return ref[]
end


end #module