using Documenter, Flight

makedocs(;
    authors="Miguel Alonso <miguel883@gmail.com>",
    sitename="Flight.jl",
    doctest = false,
    remotes = nothing,
    pages = [
        "About Flight.jl" => "index.md",
        "Showcase" => ["examples/ex01/ex01.md", "examples/ex02/ex02.md"],
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = true,
        repolink="http://localhost:8000",
    )
)

# deploydocs(; repo = "https://github.com/e271828e/Flight.jl.git")
