using Documenter
using GravityToolsLight

makedocs(
    sitename = "GravityToolsLight Documentation",
    modules = [GravityToolsLight],
    build   = "build",
    pages = [
        "Home" => "index.md",
        "Tempo Versions" => "abstract_tempo.md",
        # Add more pages as needed
    ]
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    root = joinpath(dirname(pathof(GravityToolsLight)), "..", "docs"),
    repo = "github.com/AlexBatrakov/GravityToolsLight.jl",
    target   = "dev",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages"
)
 