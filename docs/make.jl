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
    ],
    build_dir = "docs/build"
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    root = joinpath(dirname(pathof(GravityToolsLight)), "..", "docs"),
    repo = "github.com/AlexBatrakov/GravityToolsLight.jl",
    deploy_config = Documenter.GitHubActions(),
    target = "site", 
    branch = "gh-pages",
    build_dir = "docs/build",
)
 