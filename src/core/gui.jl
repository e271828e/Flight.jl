module GUI

using UnPack
using Reexport

#these are needed by any module extending GUI draw methods
@reexport using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic
@reexport using Printf
using ImGuiGLFWBackend
using ImGuiGLFWBackend.LibGLFW #defines GLFWwindow
using ImGuiGLFWBackend.LibCImGui
using ImGuiOpenGLBackend
using ImGuiOpenGLBackend.ModernGL

export CImGuiStyle, Renderer
export dynamic_button, display_bar, safe_slider, safe_input, @running_plot

################################################################################
############################# Renderer####################################

@enum CImGuiStyle begin
    classic = 0
    dark = 1
    light = 2
end

#refresh: number of display updates per frame render:
#T_render = T_display * refresh (where typically T_display = 16.67ms).
#refresh = 1 syncs the render frame rate to the display rate (vsync)
#refresh = 0 uncaps the render frame rate (WARNING: an independent scheduling
#mechanism should be used)

mutable struct Renderer
    label::String
    wsize::Tuple{Int, Int}
    style::CImGuiStyle
    refresh::Integer
    _enabled::Bool
    _initialized::Bool
    _window::Ptr{ImGuiGLFWBackend.GLFWwindow}
    _window_ctx::ImGuiGLFWBackend.Context
    _opengl_ctx::ImGuiOpenGLBackend.Context
    _cimgui_ctx::Ptr{CImGui.LibCImGui.ImGuiContext}

    function Renderer(; label = "Renderer", wsize = (1280, 720),
                        style = dark, refresh = 1)
        _enabled = true
        _initialized = false
        new(label, wsize, style, refresh, _enabled, _initialized)
    end

end

Base.propertynames(::Renderer) = (:label, :wsize, :style, :refresh)

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

    @unpack label, wsize, style, refresh, _enabled = renderer

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

    # create window
    _window = glfwCreateWindow(wsize[1], wsize[2], label, C_NULL, C_NULL)
    @assert _window != C_NULL
    glfwMakeContextCurrent(_window)
    glfwSwapInterval(refresh)

    # create OpenGL and GLFW context
    _window_ctx = ImGuiGLFWBackend.create_context(_window)
    _opengl_ctx = ImGuiOpenGLBackend.create_context()

    # setup Dear ImGui context
    _cimgui_ctx = CImGui.CreateContext()

    #comment when not using docking and multiviewports
    #enable docking and multi-viewport
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable

    # setup Dear ImGui style
    style === classic && CImGui.StyleColorsClassic()
    style === dark && CImGui.StyleColorsDark()
    style === light && CImGui.StyleColorsLight()

    #comment when not using docking and multiviewports
    # When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(style.Colors, CImGui.ImGuiCol_WindowBg, ImVec4(col.x, col.y, col.z, 1.0f0))
    end

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
        @assert renderer.refresh > 0 "The standalone run() must not be called "*
        "an unsynced Renderer (refresh = 0). Use scheduled calls to render() instead."
        #this is because a Renderer with refresh=0 does not wait for monitor refresh
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
draw(args...) = nothing

#generic mutating draw function, to be extended by users
draw!(args...) = nothing

#must be used within a CImGui.Begin() / CImGui.End() context
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
################################ Macros ########################################

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
function dynamic_button(label::String, hue::AbstractFloat)
    CImGui.PushStyleColor(CImGui.ImGuiCol_Button, CImGui.HSV(hue, 0.6, 0.6))
    CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonHovered, CImGui.HSV(hue, 0.7, 0.7))
    CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonActive, CImGui.HSV(hue, 0.8, 0.8))
    CImGui.Button((label))
    is_pressed = CImGui.IsItemActive()
    CImGui.PopStyleColor(3)
    return is_pressed
end

function display_bar(label::String, source::Real, lower_bound::Real, upper_bound::Real, size_arg = (0, 0))
    CImGui.Text(label)
    CImGui.SameLine()
    CImGui.ProgressBar((source - lower_bound)/(upper_bound - lower_bound), size_arg, "$source")
end

function safe_slider(label::String, source::AbstractFloat, lower_bound::Real, upper_bound::Real, display_format::String)
    ref = Ref(Cfloat(source))
    CImGui.Text(label)
    CImGui.SameLine()
    CImGui.SliderFloat("##"*(label), ref, lower_bound, upper_bound, display_format)
    CImGui.SameLine()
    show_help_marker("Ctrl+Click for keyboard input")
    return ref[]
end

function safe_input(label::String, source::AbstractFloat, step::Real, fast_step::Real, display_format::String)
    ref = Ref(Cdouble(source))
    CImGui.Text(label)
    CImGui.SameLine()
    CImGui.InputDouble("##"*(label), ref, step, fast_step, display_format)
    return ref[]
end

#inactive while not enabled; overwrites target while enabled
macro enabled_slider(label, target, lower_bound, upper_bound, default)
    enable = gensym(:enable)
    value = gensym(:value)
    return esc(quote
        $enable = @cstatic check=false @c CImGui.Checkbox($label, &check)
        CImGui.SameLine()
        $value = @cstatic f=Cfloat($default) @c CImGui.SliderFloat("##"*($label), &f, $lower_bound, $upper_bound)
        $enable && ($target = $value)
    end)
end

macro running_plot(label, source, lower_bound, upper_bound, initial_value, window_height)
    values = gensym(:values)
    offset = gensym(:offset)
    return esc(quote
        @cstatic $values=fill(Cfloat($initial_value),90) $offset=Cint(0) begin
            $values[$offset+1] = $source
            $offset = ($offset+1) % length($values)
            CImGui.PlotLines(string($source |> Float32), $values, length($values), $offset,
                             $label, $lower_bound, $upper_bound, (0, $window_height))
        end
    end)
end


end #module