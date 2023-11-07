using Documenter
using GravityToolsLight

makedocs(
    sitename = "GravityToolsLight Documentation",
    format = Documenter.HTML(),
    modules = [GravityToolsLight],
    pages = [
        "Home" => "index.md",
        "Tempo Versions" => "abstract_tempo.md",
        # Add more pages as needed
    ]
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    repo = "github.com/AlexBatrakov/GravityToolsLight.jl.git",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages"
)
 