using Documenter
using PenguinAnalysis

makedocs(
    modules = [PenguinAnalysis],
    authors = "PenguinxCutCell contributors",
    sitename = "PenguinAnalysis.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/PenguinAnalysis.jl",
        repolink = "https://github.com/PenguinxCutCell/PenguinAnalysis.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Examples" => "examples.md",
        "Convergence" => "convergence.md",
    ],
    pagesonly = true,
    warnonly = false,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/PenguinAnalysis.jl",
        push_preview = true,
    )
end
