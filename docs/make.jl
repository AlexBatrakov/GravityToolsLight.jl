using Documenter
using GravityToolsLight

- name: Install Python dependencies
  run: |
    julia --project=docs/ -e 'using Pkg; Pkg.add("Conda"); using Conda; Conda.add("matplotlib")'


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
    branch = "gh-pages",
    deploy_config = Documenter.GitHubActions(),
    push_preview = true, # set to false when ready for production
)
 