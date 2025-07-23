using Documenter, Flight

makedocs(;
    authors="Miguel Alonso <miguel883@gmail.com>",
    sitename="Flight.jl",
    # repo="https://github.com/e271828e/Flight.jl.git",
    doctest = false,
    # remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Showcase" => ["examples/ex01/ex01.md", "examples/ex02/ex02.md"],
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = true,
        # repolink="https://github.com/e271828e/Flight.jl.git",
    )
)

deploydocs(; repo = "https://github.com/e271828e/Flight.jl.git")
