using CImGui
using CImGui.CSyntax
using CImGui.CSyntax.CStatic
using CImGui.GLFWBackend
using CImGui.OpenGLBackend
using CImGui.GLFWBackend.GLFW
using CImGui.OpenGLBackend.ModernGL
using Printf

Base.@kwdef struct DummySim
    t::Ref{Float64} = Ref(0.0)
    y::Ref{Float64} = Ref(0.0)
end

function step!(sim::DummySim, dt::Float64)
    f = 0.2
    sim.t[] += dt
    sim.y[] = sin(2π*f*sim.t[])
end

reset!(sim::DummySim) = (sim.t[] = 0)

function get_τ_source()
    let wall_time_ref = time()
        ()-> time() - wall_time_ref
    end
end

function sim_loop(sim::DummySim, channel::Channel{Float64},
                τ::Function; dt::Real = 0.02)

    reset!(sim)
    steps = 0

    println("Sim loop starting at thread $(Threads.threadid())...")
    while true

        #problem: at the time we enter this loop iteration, render_task might
        #not be done yet, but when we get to put!, it is, so there is no one to
        #take! it, and this task gets permanently blocked. therefore, we can't
        #rely on while !istaskdone(render_task) to know when we have to exit
        #this loop. instead, we bind the channel to the render task, and catch
        #the exception that will . if the channel is open, the render task has
        #not terminated yet (however, there is a small chance that the channel )

        # τ_before = τ()
        println("Took $steps steps, putting y(t = $(sim.t[])) = $(sim.y[]) at τ = $(τ())")
        try
            put!(channel, sim.y[]) #will block while the buffer is full
        catch err
            if err isa InvalidStateException
                println("Channel closed, sim loop exiting...")
                break
            else
                rethrow()
            end
        end
        # Δτ_put = τ() - τ_before
        # println("Time blocked at put! = $Δτ_put")

        steps = 0
        while sim.t[] <= τ()
            step!(sim, dt)
            steps += 1
        end
    end

end


# must use Threads.@spawn, not @spawn, which is a deprecated version of @spawnat
function start_loops_threaded(sim::DummySim = DummySim())

    @show Threads.nthreads()
    τ = get_τ_source()
    channel = Channel{Float64}()

    Threads.@spawn sim_loop(sim, channel, τ)

    wait(Threads.@spawn render_loop(channel, τ))
    close(channel)

    return render_task, sim_task, channel

end

function start_loops(sim::DummySim = DummySim())

    τ = get_τ_source()
    channel = Channel{Float64}()

    render_task = @task render_loop(channel, τ)
    bind(channel, render_task)

    sim_task = @task sim_loop(sim, channel, τ)

    schedule(sim_task)
    schedule(render_task)

    wait(render_task)

    return render_task, sim_task, channel

end

#parts of this are general, parts are output value-specific. we should extract
#the value-specific ones into a method draw(::System{<:ComponentY}). Problem is
#the type of ac.y is a horrendous NamedTuple. maybe we should not pass
#everything to these functions. maybe we should define an immutable struct
#C172R.DisplayOutput, constructable from y, which packs all the variables of
#interest into a single struct. however, we could also define a method
#draw(::System{<:Cessna172R}, y) and use it for dispatch. then we simply need
#to dispatch as draw(sim.sys.aircraft, sim.y.aircraft). this way, we don't have
#to extract the output type of the aircraft
function render_loop(channel::Channel{Float64}, τ::Function)

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
    GLFW.SwapInterval(1)  # enable vsync

    # setup Dear ImGui context
    ctx = CImGui.CreateContext()

    # setup Dear ImGui style
    CImGui.StyleColorsDark()
    # CImGui.StyleColorsClassic()
    # CImGui.StyleColorsLight()

    # load Fonts
    # - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use `CImGui.PushFont/PopFont` to select them.
    # - `CImGui.AddFontFromFileTTF` will return the `Ptr{ImFont}` so you can store it if you need to select the font among multiple.
    # - If the file cannot be loaded, the function will return C_NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    # - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling `CImGui.Build()`/`GetTexDataAsXXXX()``, which `ImGui_ImplXXXX_NewFrame` below will call.
    # - Read 'fonts/README.txt' for more instructions and details.
    fonts_dir = joinpath(@__DIR__, "..", "fonts")
    fonts = CImGui.GetIO().Fonts
    # default_font = CImGui.AddFontDefault(fonts)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Cousine-Regular.ttf"), 15)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "DroidSans.ttf"), 16)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Karla-Regular.ttf"), 10)
    # CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "ProggyTiny.ttf"), 10)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "Roboto-Medium.ttf"), 16)
    # @assert default_font != C_NULL

    # setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true)
    ImGui_ImplOpenGL3_Init(glsl_version)

    try
        show_demo_window = true
        show_another_window = false
        clear_color = Cfloat[0.45, 0.55, 0.60, 1.00]


        println("Render loop starting at thread $(Threads.threadid())...")
        while !GLFW.WindowShouldClose(window)


            # τ_start = τ()
            value = take!(channel)
            println("Took value $value at τ = $(τ())")

            GLFW.PollEvents()
            # start the Dear ImGui frame
            ImGui_ImplOpenGL3_NewFrame()
            ImGui_ImplGlfw_NewFrame()
            CImGui.NewFrame()

            # # show the big demo window
            # show_demo_window && @c CImGui.ShowDemoWindow(&show_demo_window)

            # show a simple window that we create ourselves.
            # we use a Begin/End pair to created a named window.
            @cstatic f=Cfloat(0.0) counter=Cint(0) begin
                CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
                CImGui.Text("This is some useful text.")  # display some text
                @c CImGui.Checkbox("Demo Window", &show_demo_window)  # edit bools storing our window open/close state
                @c CImGui.Checkbox("Another Window", &show_another_window)

                @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a slider from 0 to 1
                CImGui.ColorEdit3("clear color", clear_color)  # edit 3 floats representing a color
                CImGui.Button("Button") && (counter += 1)

                CImGui.SameLine()
                CImGui.Text("counter = $counter")
                CImGui.Text("value = $value")
                CImGui.Text(@sprintf("Application average %.3f ms/frame (%.1f FPS)", 1000 / CImGui.GetIO().Framerate, CImGui.GetIO().Framerate))

                CImGui.End()
            end

            # show another simple window.
            if show_another_window
                @c CImGui.Begin("Another Window", &show_another_window)  # pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
                CImGui.Text("Hello from another window!")
                CImGui.Button("Close Me") && (show_another_window = false;)
                CImGui.End()
            end

            # rendering
            CImGui.Render()
            GLFW.MakeContextCurrent(window)
            display_w, display_h = GLFW.GetFramebufferSize(window)
            glViewport(0, 0, display_w, display_h)
            glClearColor(clear_color...)
            glClear(GL_COLOR_BUFFER_BIT)
            ImGui_ImplOpenGL3_RenderDrawData(CImGui.GetDrawData())

            GLFW.MakeContextCurrent(window)
            GLFW.SwapBuffers(window)

            # Δτ = τ() - τ_start
            # println("Frame render took Δτ = $Δτ")
        end
    catch e
        @error "Error in renderloop!" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
    finally
        ImGui_ImplOpenGL3_Shutdown()
        ImGui_ImplGlfw_Shutdown()
        CImGui.DestroyContext(ctx)
        GLFW.DestroyWindow(window)

    end

    println("Render loop exiting...")

    #leave the channel empty before exiting, in case the sim_loop has put! a
    #value since the beginning of the last render loop iteration
    isready(channel) ? take!(channel) : nothing
    return nothing

end
