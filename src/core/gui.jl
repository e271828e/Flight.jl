module GUI

using UnPack
using Reexport
using StaticArrays
using Logging

@reexport using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic
@reexport using Printf

using GLFW
using ModernGL
using CImGui.lib

using ..IODevices

function get_glsl_version(gl_version)
    gl2glsl = Dict(v"2.0" => 110, v"2.1" => 120, v"3.0" => 130, v"3.1" => 140, v"3.2" => 150)
    if gl_version < v"3.3"
        gl2glsl[gl_version]
    else
        gl_version.major * 100 + gl_version.minor * 10
    end
end

#if sync > 0:
#T_render = T_display * sync (where typically T_display = 16.67ms).
#sync = 1 syncs the refresh rate to the display rate (vsync)
#if sync = 0:
#uncaps the refresh rate (to be used only with scheduled calls to render())

mutable struct Renderer{T} <: IODevice{T}
    label::String
    monitor_pref::UInt8 #preferred monitor when multiple monitors available
    font_size::UInt8 #will be scaled by the display's content scale
    sync::UInt8 #number of display updates per frame render
    f_draw::Function #GUI function to be called
    _initialized::Bool
    _ctx::Ptr{CImGui.lib.ImGuiContext}
    _window::GLFW.Window

    function Renderer(; label = "Renderer", monitor_pref = 2, font_size = 16,
        sync = 1, f_draw = ()->nothing)
        _initialized = false
        new{Nothing}(label, monitor_pref, font_size, sync, f_draw, _initialized)
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

end

function gui_test()
    # setup Dear ImGui context
    ctx = CImGui.CreateContext()
    @show typeof(ctx)

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

    monitor_pref = 1 ##################################### UPDATE THIS
    available_monitors = GLFW.GetMonitors()
    monitor_pref = min(monitor_pref, length(available_monitors))
    monitor = available_monitors[monitor_pref]
    vmode = GLFW.GetVideoMode(monitor)

    window_size=(vmode.width, vmode.height)
    window_title="CImGui"

    # x_scale, y_scale = GLFW.GetMonitorContentScale(monitor)
    # font_scaling = max(x_scale, y_scale)
    # font_size = 12 * font_scaling

    font_size = 12 * vmode.height / 1080
    fonts_dir = joinpath(@__DIR__, "gui", "fonts")
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    @assert (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Sans Linear-Regular.ttf"), font_size) != C_NULL)
    #@assert (CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Mono Linear-Regular.ttf"), font_size) != C_NULL)

    opengl_version = v"3.2" #Only versions from 3.2 to 4.1 currently work on MacOS

    # Configure GLFW
    glsl_version = get_glsl_version(opengl_version)
    GLFW.WindowHint(GLFW.VISIBLE, true)
    GLFW.WindowHint(GLFW.DECORATED, true)
    GLFW.WindowHint(GLFW.FOCUSED, true)
    GLFW.WindowHint(GLFW.MAXIMIZED, false)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, opengl_version.major)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, opengl_version.minor)

    if Sys.isapple()
        GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE)
        GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, ModernGL.GL_TRUE)
    end

    # Create window
    global _window = GLFW.CreateWindow(window_size[1], window_size[2], window_title)
    window = _window
    @assert window != C_NULL

    @show typeof(window)

    # _, y_pos = GLFW.GetWindowPos(_window)
    # GLFW.SetWindowPos(_window, 0, y_pos) #no effect on borderless window

    GLFW.MakeContextCurrent(window)

    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############
    GLFW.SwapInterval(1)  # enable vsync
    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############
    ############# FIX THIS ############

    # Setup Platform/Renderer bindings
    CImGui.lib.ImGui_ImplGlfw_InitForOpenGL(Ptr{CImGui.lib.GLFWwindow}(window.handle), true)
    CImGui.lib.ImGui_ImplOpenGL3_Init("#version $(glsl_version)")

    try
        while !GLFW.WindowShouldClose(window)
            GLFW.PollEvents()

            # Start the Dear ImGui frame
            CImGui.lib.ImGui_ImplOpenGL3_NewFrame()
            CImGui.lib.ImGui_ImplGlfw_NewFrame()
            CImGui.NewFrame()

            #DRAW FUNCTION GOES HERE
            #DRAW FUNCTION GOES HERE
            #DRAW FUNCTION GOES HERE
            @cstatic f=Cfloat(0.0) counter=Cint(0) begin
                CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
                CImGui.Text("This is some useful text.")  # display some text

                @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a slider from 0 to 1
                CImGui.Button("Button") && (counter += 1)

                CImGui.SameLine()
                CImGui.Text("counter = $counter")
                CImGui.Text(@sprintf("Application average %.3f ms/frame (%.1f FPS)", 1000 / unsafe_load(CImGui.GetIO().Framerate), unsafe_load(CImGui.GetIO().Framerate)))

                CImGui.End()
            end
            #DRAW FUNCTION GOES HERE
            #DRAW FUNCTION GOES HERE
            #DRAW FUNCTION GOES HERE
            #DRAW FUNCTION GOES HERE

            # Rendering
            CImGui.Render()
            GLFW.MakeContextCurrent(window)

            display_w, display_h = GLFW.GetFramebufferSize(window)

            ModernGL.glViewport(0, 0, display_w, display_h)
            ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT)
            CImGui.lib.ImGui_ImplOpenGL3_RenderDrawData(Ptr{Cint}(CImGui.GetDrawData()))

            GLFW.MakeContextCurrent(window)
            GLFW.SwapBuffers(window)

            if (unsafe_load(CImGui.lib.igGetIO().ConfigFlags) & CImGui.lib.ImGuiConfigFlags_ViewportsEnable) == CImGui.lib.ImGuiConfigFlags_ViewportsEnable
                backup_current_context = GLFW.GetCurrentContext()
                CImGui.lib.igUpdatePlatformWindows()
                CImGui.lib.igRenderPlatformWindowsDefault(C_NULL, C_NULL)
                GLFW.MakeContextCurrent(backup_current_context)
            end

        end
    catch e
        @error "Error in CImGui $(CImGui._backend[]) renderloop!" exception=(e, catch_backtrace())
    finally
        CImGui.lib.ImGui_ImplOpenGL3_Shutdown()
        CImGui.lib.ImGui_ImplGlfw_Shutdown()
        CImGui.DestroyContext(ctx)
        GLFW.DestroyWindow(window)
    end

    # CImGui.render(ctx; engine, clear_color=Ref(clear_color)) do
    #     if isnothing(image_id)
    #         image_id = CImGui.create_image_texture(img_width, img_height)
    #     end

    #     demo_open && @c ShowJuliaDemoWindow(&demo_open)

    #     # show image example
    #     if CImGui.Begin("Image Demo")
    #         image = rand(ModernGL.GLubyte, 4, img_width, img_height)
    #         CImGui.update_image_texture(image_id, image, img_width, img_height)
    #         CImGui.Image(Ptr{Cvoid}(image_id), CImGui.ImVec2(img_width, img_height))
    #         CImGui.End()
    #     end
    # end
end


end