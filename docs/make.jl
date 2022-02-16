using Documenter
# using Dojo  # your package name here

 makedocs(
     sitename = "Dojo",  # your package name here
     format = Documenter.HTML(prettyurls = false),  # optional
     pages = [
         "Introduction" => "index.md",
         "Examples" => "examples.md",
         "Notations" => "notations.md"
     ]
 )

 # Documenter can also automatically deploy documentation to gh-pages.
 # See "Hosting Documentation" and deploydocs() in the Documenter manual
 # for more information.
 deploydocs(
     repo = "github.com/dojo-sim/Dojo.jl.git",
 )