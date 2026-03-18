module GUI

using Reexport
using StaticArrays
using Logging

@reexport using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic
@reexport using Printf

using GLFW
using ModernGL
using CImGui.lib

using ..IODevices

export CImGuiStyle, Renderer
export mode_button, dynamic_button, toggle_switch, display_bar, safe_slider, safe_input
export HSV_amber, HSV_gray, HSV_green, HSV_red


################################################################################
############################### Renderer #######################################

#sync > 0: syncs the renderer's refresh rate to the display rate. more precisely
#T_render = T_display * sync (where often T_display = 16.67ms).

#sync = 0: uncaps the refresh rate (should only be used with independently
#scheduled calls to render())

mutable struct Renderer <: IODevice
    label::String #window label
    monitor_pref::UInt8 #preferred monitor when multiple monitors available
    font_size::UInt8 #will be scaled by the display's content scale
    sync::UInt8 #see above
    f_draw::Function #GUI IO function to be called
    _initialized::Bool
    _ctx::Ptr{CImGui.lib.ImGuiContext}
    _window::GLFW.Window

    function Renderer(; label = "Renderer", monitor_pref = 2, font_size = 16,
        sync = 1, f_draw = ()->nothing)
        _initialized = false
        new(label, monitor_pref, font_size, sync, f_draw, _initialized)
    end

end

Base.propertynames(::Renderer) = (:label, :monitor_pref, :font_size, :sync, :f_draw)

function Base.setproperty!(renderer::Renderer, name::Symbol, value)
    if name ∈ propertynames(renderer)
        if renderer._initialized
            @error("Cannot set property $name for an initialized Renderer, ",
            "call shutdown! first")
        else
            setfield!(renderer, name, value)
        end
    else
        @error("Unsupported property: $name")
    end
end


function IODevices.init!(renderer::Renderer)

    (; label, monitor_pref, font_size, sync) = renderer

    # setup Dear ImGui context
    _ctx = CImGui.CreateContext()

    # enable docking and multi-viewport
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable

    # When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(style.Colors, CImGui.ImGuiCol_WindowBg, ImVec4(col.x, col.y, col.z, 1.0f0))
    end

    # setup Dear ImGui style
    CImGui.StyleColorsDark()
    # CImGui.StyleColorsClassic()
    # CImGui.StyleColorsLight()

    available_monitors = GLFW.GetMonitors()
    monitor_pref = min(monitor_pref, length(available_monitors))
    monitor = available_monitors[monitor_pref]
    vmode = GLFW.GetVideoMode(monitor)
    window_size=(vmode.width, vmode.height)

    font_size = 12 * vmode.height / 1080
    fonts_dir = joinpath(@__DIR__, "gui", "fonts")
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    @assert (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Sans Linear-Regular.ttf"), font_size) != C_NULL)
    #@assert (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Mono Linear-Regular.ttf"), font_size) != C_NULL)

    #set OpenGL version and its corresponding GLSL version
    gl_version = v"3.2" #Only versions from 3.2 to 4.1 currently work on MacOS
    gl2glsl = Dict(v"2.0" => 110, v"2.1" => 120, v"3.0" => 130, v"3.1" => 140, v"3.2" => 150)
    glsl_version = (gl_version < v"3.3" ? gl2glsl[gl_version] : 100gl_version.major + 10gl_version.minor)

    GLFW.WindowHint(GLFW.VISIBLE, true)
    GLFW.WindowHint(GLFW.DECORATED, true)
    GLFW.WindowHint(GLFW.FOCUSED, true)
    GLFW.WindowHint(GLFW.MAXIMIZED, false)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, gl_version.major)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, gl_version.minor)

    if Sys.isapple()
        GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE)
        GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, ModernGL.GL_TRUE)
    end

    _window = GLFW.CreateWindow(window_size[1], window_size[2], label)
    @assert _window != C_NULL

    GLFW.MakeContextCurrent(_window)
    GLFW.SwapInterval(sync)  # enable vsync

    #setup Platform/Renderer bindings
    CImGui.lib.ImGui_ImplGlfw_InitForOpenGL(Ptr{CImGui.lib.GLFWwindow}(_window.handle), true)
    CImGui.lib.ImGui_ImplOpenGL3_Init("#version $(glsl_version)")

    setfield!(renderer, :_ctx, _ctx)
    setfield!(renderer, :_window, _window)
    setfield!(renderer, :_initialized, true)

    return nothing

end

function IODevices.shutdown!(renderer::Renderer)

    @assert renderer._initialized "Cannot shutdown an uninitialized renderer"

    CImGui.lib.ImGui_ImplOpenGL3_Shutdown()
    CImGui.lib.ImGui_ImplGlfw_Shutdown()
    CImGui.DestroyContext(renderer._ctx)
    GLFW.DestroyWindow(renderer._window)
    setfield!(renderer, :_initialized, false)

    return nothing

end

function IODevices.should_close(renderer::Renderer)
    renderer._initialized ? Bool(GLFW.WindowShouldClose(renderer._window)) : false
end

function render!(renderer::Renderer)

    (; f_draw, _initialized, _window) = renderer

    @assert _initialized "Renderer not initialized, call init! before update!"

    GLFW.PollEvents()

    # Start the Dear ImGui frame
    CImGui.lib.ImGui_ImplOpenGL3_NewFrame()
    CImGui.lib.ImGui_ImplGlfw_NewFrame()
    CImGui.NewFrame()

    #call user-defined IO function
    f_draw()

    # Rendering
    CImGui.Render()
    GLFW.MakeContextCurrent(_window)

    display_w, display_h = GLFW.GetFramebufferSize(_window)

    ModernGL.glViewport(0, 0, display_w, display_h)
    ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT)
    CImGui.lib.ImGui_ImplOpenGL3_RenderDrawData(Ptr{Cint}(CImGui.GetDrawData()))

    GLFW.MakeContextCurrent(_window)
    GLFW.SwapBuffers(_window)

    if (unsafe_load(CImGui.GetIO().ConfigFlags) & CImGui.lib.ImGuiConfigFlags_ViewportsEnable) == CImGui.lib.ImGuiConfigFlags_ViewportsEnable
        backup_current_context = GLFW.GetCurrentContext()
        CImGui.lib.igUpdatePlatformWindows()
        CImGui.lib.igRenderPlatformWindowsDefault(C_NULL, C_NULL)
        GLFW.MakeContextCurrent(backup_current_context)
    end

end

function render_loop(renderer::Renderer, timeout::Real = Inf64)

    try
        renderer._initialized || IODevices.init!(renderer)
    catch e
        @warn "Error while initializing Renderer" exception=e
    end

    renderer._initialized || return

    try
        @assert renderer.sync > 0 "The standalone render_loop() must not be called "*
        "for an unsynced Renderer (sync = 0). Use scheduled calls to update!() instead."
        #a Renderer with sync=0 does not wait for monitor sync when
        #glfwSwapBuffers is called within update!(). this means the GUI frame
        #rate is effectively uncapped. this is generally not good, so the calls
        #to update! must be frequency-limited by some other scheduling means

        t0 = time()
        while !IODevices.should_close(renderer) && (time() - t0 < timeout)
            render!(renderer)
        end
    catch e
        @warn "Error while updating window" exception=e
        # Base.show_backtrace(stderr, catch_backtrace())
    finally
        IODevices.shutdown!(renderer)
    end

end

#generic mutating draw function, to be extended by users. if not specialized,
#falls back to non-mutating version
draw!(args...; kwargs...) = draw(args...; kwargs...)

#generic non-mutating frame draw function, to be extended by users
draw(args...; kwargs...) = nothing

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


################################################################################
########################### Custom Widgets #####################################

const HSV_gray = (0.0, 0.0, 0.3)
const HSV_amber = (0.13, 0.6, 0.6)
const HSV_green = (0.4, 0.6, 0.6)
const HSV_red = (0.0, 0.7, 0.7)

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
                        idle_HSV::NTuple{3,Real};
                        ΔSV_hover::Real = 0.1,
                        ΔSV_push::Real = 0.2)

    hover_HSV = (idle_HSV[1], idle_HSV[2] + ΔSV_hover, idle_HSV[3] + ΔSV_hover)
    push_HSV = (idle_HSV[1], idle_HSV[2] + ΔSV_push, idle_HSV[3] + ΔSV_push)
    dynamic_button(label, idle_HSV, hover_HSV, push_HSV)

end

function mode_button(label::String,
                    button_mode::T, #mode enabled by this button
                    requested_mode::T, #currently requested mode
                    active_mode::T; #currently active mode
                    HSV_none::NTuple{3,Real} = HSV_gray,
                    HSV_requested::NTuple{3,Real} = HSV_amber,
                    HSV_active::NTuple{3,Real} = HSV_green, kwargs...) where {T}

    if active_mode === button_mode #mode enabled by this button is active
        HSV_button = HSV_active
    elseif requested_mode === button_mode #mode enabled by this button is requested but not active
        HSV_button = HSV_requested
    else #mode enabled by this button is neither active nor requested
        HSV_button = HSV_none
    end

    dynamic_button(label, HSV_button, kwargs...)

end

function display_bar(label::String, current_value::Real, lower_bound::Real, upper_bound::Real, size_arg = (0, 0))
    CImGui.Text(label)
    CImGui.SameLine()
    CImGui.ProgressBar((current_value - lower_bound)/(upper_bound - lower_bound), size_arg, "$current_value")
end

#the string after ## is not shown, but is part of the widget's ID, which must be
#unique to avoid conflicts. an alternative solution is to use PushID and PopID.
#See DearImGui FAQ
function safe_slider(label::String, current_value::AbstractFloat, args...)
    ref = Ref(Cfloat(current_value))
    CImGui.SliderFloat(label, ref, args...)
    # CImGui.SameLine()
    # show_help_marker("Ctrl+Click for keyboard input")
    return ref[]
end

#ref is a stack-allocated variable. we return the value it points to, not the
#Ref itself, so no leftover memory
function safe_input(label::String, current_value::AbstractFloat, args...)
    ref = Ref(Cdouble(current_value))
    CImGui.InputDouble(label, ref, args...)
    return ref[]
end

end #module