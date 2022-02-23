## docs/make.jl --- to make the docs

### got this line from DynamicalSystems.jl
CI = get(ENV, "CI", nothing) === true || !isnothing(get(ENV, "GITHUB_TOKEN", nothing))

using Documenter
using ProtStructs

makedocs(
    sitename="ProtStructs",
    authors="Bon Leif Amalla and contributors",
    doctest=false,
    format = Documenter.HTML(prettyurls=CI),
    pages= [
        "Introduction" => "index.md",
        # "Getting Started" => "start.md",
        # "Hydrogen bonding model" => "hbond.md",
        # "Secondary structure" => "secstruct.md",
        # "Implementation notes" => "impnotes.md"
   ]
)

if CI
    deploydocs(
        repo = "https://github.com/bldamalla/ProtStructs.jl.git",
        target = "build",
        push_preview = false
    )
end

