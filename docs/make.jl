using Documenter
using GravityToolsLight

makedocs(
    root = joinpath(dirname(pathof(GravityToolsLight)), "..", "docs"),
    sitename = "GravityToolsLight Documentation",
    modules = [GravityToolsLight],
    pages = [
        "Home" => "index.md",
        "Tempo Versions" => "abstract_tempo.md",
        # Add more pages as needed
    ]
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    repo = "github.com/AlexBatrakov/GravityToolsLight.jl",
    target = "dev"
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages"
)
 