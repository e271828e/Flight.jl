

function test1()
    @macroexpand begin
    @cstatic f=Cfloat(0.0) begin
        CImGui.Begin("Hello, world!")  # create a window called "Hello, world!" and append into it.
        CImGui.Text("This is some useful text.")  # display some text

        @c CImGui.SliderFloat("float", &f, 0, 1)  # edit 1 float using a #         slider from 0 to 1

        CImGui.End()
    end

end

end