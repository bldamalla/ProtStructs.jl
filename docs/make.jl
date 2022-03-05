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
    modules = [ProtStructs,],
    pages= [
        "Introduction" => "index.md",
        "Brief walkthrough" => "start.md",
        "Models" => [
            "Structure Frames" => "structs/sfs.md",
            # "HBond Dictionary" => "structs/hbdict.md"
        ],
        # "Hydrogen bonding model" => "hbond.md",
        # "Secondary structure" => "secstruct.md",
        # "Implementation notes" => "impnotes.md",
        "API Documentation" => "API.md"
   ]
)

if CI
    deploydocs(
        repo = "https://github.com/bldamalla/ProtStructs.jl.git",
        target = "build",
        push_preview = false
    )
end

