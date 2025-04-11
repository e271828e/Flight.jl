using Documenter, Flight

makedocs(;
    authors="Miguel Alonso <miguel883@gmail.com>",
    sitename="Flight.jl",
    doctest = false,
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = true,
        repolink="http://localhost:8000",
    )
)