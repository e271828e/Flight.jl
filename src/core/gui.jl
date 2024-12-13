module GUI

using Printf
using CImGui
using CImGui.lib
import CImGui.CSyntax: @c, @cstatic

# Load deps for the GLFW/OpenGL backend
import GLFW
import ModernGL

export gui_test

function get_glsl_version(gl_version)
    gl2glsl = Dict(v"2.0" => 110, v"2.1" => 120, v"3.0" => 130, v"3.1" => 140, v"3.2" => 150)
    if gl_version < v"3.3"
        gl2glsl[gl_version]
    else
        gl_version.major * 100 + gl_version.minor * 10
    end
end

function gui_test()
    # setup Dear ImGui context
    ctx = CImGui.CreateContext()

    # enable docking and multi-viewport
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable

    # setup Dear ImGui style
    CImGui.StyleColorsDark()
    # CImGui.StyleColorsClassic()
    # CImGui.StyleColorsLight()

    # When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(style.Colors, CImGui.ImGuiCol_WindowBg, ImVec4(col.x, col.y, col.z, 1.0f0))
    end

    # load Fonts
    # - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use `CImGui.PushFont/PopFont` to select them.
    # - `CImGui.AddFontFromFileTTF` will return the `Ptr{ImFont}` so you can store it if you need to select the font among multiple.
    # - If the file cannot be loaded, the function will return C_NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    # - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling `CImGui.Build()`/`GetTexDataAsXXXX()``, which `ImGui_ImplXXXX_NewFrame` below will call.
    # - Read 'fonts/README.txt' for more instructions and details.
    fonts_dir = joinpath(@__DIR__, "..", "fonts")
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    # default_font = CImGui.AddFontDefault(fonts)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Cousine-Regular.ttf"), 15)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "DroidSans.ttf"), 16)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Karla-Regular.ttf"), 10)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "ProggyTiny.ttf"), 10)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Roboto-Medium.ttf"), 16)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Mono Casual-Regular.ttf"), 16)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Mono Linear-Regular.ttf"), 16)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Sans Casual-Regular.ttf"), 16)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Recursive Sans Linear-Regular.ttf"), 16)
    # @assert default_font != C_NULL

    ############### Backend-specific stuff starts here #########
    ############### Backend-specific stuff starts here #########
    ############### Backend-specific stuff starts here #########
    ############### Backend-specific stuff starts here #########
    opengl_version = v"3.2" #lowest that will work on MacOS, highest is currently 4.1
    window_size=(1280, 720)
    window_title="CImGui"

    # Configure GLFW
    glsl_version = get_glsl_version(opengl_version)
    GLFW.WindowHint(GLFW.VISIBLE, true)
    GLFW.WindowHint(GLFW.DECORATED, true)
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