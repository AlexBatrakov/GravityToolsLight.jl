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
    ],
    checkdocs = :none
)

# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    repo = "https://github.com/AlexBatrakov/GravityToolsLight.jl.git",
    deploy_config = Documenter.GitHubActions(),
    target_branch="gh-pages",
    push_preview = true, # set to false when ready for production
)
 