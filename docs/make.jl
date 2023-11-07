using Documenter
using GravityToolsLight

makedocs(
    sitename = "GravityToolsLight Documentation",
    format = Documenter.HTML(),
    modules = [GravityToolsLight],
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Installation" => "installation.md",
        "Quick Start" => "quickstart.md",
        
        "GeneralTempoFramework" => [
            "Overview" => "frameworks/tempo/overview.md",
            "Getting Started" => "frameworks/tempo/getting_started.md",
            "Tempo Versions" => "frameworks/tempo/abstract_tempo.md",
            "API Reference" => "frameworks/tempo/api_reference.md",
            "Examples" => "frameworks/tempo/examples.md"
        ],
        
        "GeneralPKFramework" => [
            "Overview" => "frameworks/pk/overview.md",
            "Getting Started" => "frameworks/pk/getting_started.md",
            "API Reference" => "frameworks/pk/api_reference.md",
            "Examples" => "frameworks/pk/examples.md"
        ],
        
        "AdaptiveRefinement2DGrid" => [
            "Overview" => "frameworks/refinement/overview.md",
            "Getting Started" => "frameworks/refinement/getting_started.md",
            "API Reference" => "frameworks/refinement/api_reference.md",
            "Examples" => "frameworks/refinement/examples.md"
        ],
        
        "Combining Frameworks" => [
            "Tempo and Refinement" => "combining/tempo_refinement.md",
            "PK and Refinement" => "combining/pk_refinement.md",
            "Tempo, PK, and Refinement" => "combining/all_three.md"
        ],
        
        "Additional Information" => [
            "FAQ" => "additional/faq.md",
            "Contributing" => "additional/contributing.md",
            "License" => "additional/license.md"
        ],
        
        "Technical Details" => [
            "Architecture" => "technical/architecture.md",
            "Performance" => "technical/performance.md",
            "Changelog" => "technical/changelog.md"
        ],
        
        "Glossary" => "glossary.md",
        "Index" => "index.md"
    ],
    # Other configurations as needed
)



# Optionally, include to push the docs to the gh-pages branch
deploydocs(
    repo = "github.com/AlexBatrakov/GravityToolsLight.jl.git",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages"
)
 