module GUI

using UnPack
using Reexport

using CImGui.GLFWBackend
using CImGui.OpenGLBackend
using CImGui.GLFWBackend.GLFW
using CImGui.OpenGLBackend.ModernGL

#these are needed by any module extending GUI draw methods
@reexport using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic
@reexport using Printf

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
#refresh = 0 uncaps the render frame rate (use judiciously!)

mutable struct Renderer
    label::String
    wsize::Tuple{Int, Int}
    style::CImGuiStyle
    refresh::Integer
    _enabled::Bool
    _initialized::Bool
    _window::GLFW.Window
    _context::Ptr{CImGui.LibCImGui.ImGuiContext}

    function Renderer(; label = "Renderer", wsize = (1280, 720),
                        style = dark, refresh = 0)
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
    _window = GLFW.CreateWindow(wsize[1], wsize[2], label)
    @assert _window != C_NULL
    GLFW.MakeContextCurrent(_window)
    GLFW.SwapInterval(refresh)

    # setup Dear ImGui context
    _context = CImGui.CreateContext()

    # setup Dear ImGui style
    style === classic && CImGui.StyleColorsClassic()
    style === dark && CImGui.StyleColorsDark()
    style === light && CImGui.StyleColorsLight()

    # setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(_window, true)
    ImGui_ImplOpenGL3_Init(glsl_version)

    setfield!(renderer, :_initialized, true)
    setfield!(renderer, :_window, _window)
    setfield!(renderer, :_context, _context)

    return nothing

end


function render(renderer::Renderer, fdraw!::Function, fdraw_args...)

    renderer._enabled || return

    @unpack _initialized, _window = renderer

    @assert _initialized "Renderer not initialized, call init! before update!"

    try
        # start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame()
        ImGui_ImplGlfw_NewFrame()
        CImGui.NewFrame()

        #draw the frame and apply user inputs to arguments
        fdraw!(fdraw_args...)

        CImGui.Render()
        GLFW.MakeContextCurrent(_window)

        display_w, display_h = GLFW.GetFramebufferSize(_window)
        glViewport(0, 0, display_w, display_h)
        glClear(GL_COLOR_BUFFER_BIT)
        ImGui_ImplOpenGL3_RenderDrawData(CImGui.GetDrawData())

        GLFW.MakeContextCurrent(_window)
        GLFW.SwapBuffers(_window)
        GLFW.PollEvents() #essential to catch window close requests

    catch e

        @error "Error while updating window" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
        shutdown!(renderer)

    end

    return nothing

end


function run(renderer::Renderer, fdraw!::Function, fdraw_args...)

    renderer._enabled || return
    renderer._initialized || init!(renderer)

    while !GLFW.WindowShouldClose(renderer._window)
        render(renderer, fdraw!, fdraw_args...)
    end

    shutdown!(renderer)

end


function should_close(renderer::Renderer)

    renderer._enabled || return false
    renderer._initialized ? GLFW.WindowShouldClose(renderer._window) : false

end


function shutdown!(renderer::Renderer)

    renderer._enabled || return
    @assert renderer._initialized "Cannot shutdown an uninitialized renderer"

    ImGui_ImplOpenGL3_Shutdown()
    ImGui_ImplGlfw_Shutdown()
    CImGui.DestroyContext(renderer._context)
    GLFW.DestroyWindow(renderer._window)
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

function draw_test()

    output = @cstatic f=Cfloat(0.0) begin
        CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
        CImGui.Text("This is some useful text.")  # display some text
        @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a slider from 0 to 1
        CImGui.End()
    end

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