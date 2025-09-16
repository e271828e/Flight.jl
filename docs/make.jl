using Documenter, Flight

makedocs(;
    authors="Miguel Alonso <miguel883@gmail.com>",
    sitename="Flight.jl",
    # repo="https://github.com/e271828e/Flight.jl.git",
    doctest = false,
    # remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/tutorial01/tutorial01.md",
            "tutorials/tutorial02/tutorial02.md",
            ],
        # "API" => "api.md",
    ],
    draft = false,
    pagesonly = false,
    format = Documenter.HTML(
        prettyurls = true,
        size_threshold = 1000*1024,
        size_threshold_warn = 500*1024,
        example_size_threshold = 100*1024,
        # repolink="https://github.com/e271828e/Flight.jl.git",
    )
)

deploydocs(; repo = "https://github.com/e271828e/Flight.jl.git")
