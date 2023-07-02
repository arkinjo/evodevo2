using Pkg
Pkg.add("PkgTemplates")

using PkgTemplates

t = Template(;
           user="arkinjo",
           authors=["Akira Kinjo"],
           plugins=[
               License(name="CC0"),
               Git(),
               GitHubActions(),
           ],
       )
