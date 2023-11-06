using Documenter
using GravityToolsLight

makedocs(
    sitename = "GravityToolsLight Documentation",
    modules = [GravityToolsLight],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Tempo Versions" => "abstract_tempo.md",
        # Add more pages as needed
    ]
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    repo = "https://github.com/AlexBatrakov/GravityToolsLight.jl.git",
    push_preview = true, # set to false when ready for production
)