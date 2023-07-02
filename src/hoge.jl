using Pkg
Pkg.add("PkgTemplates")

using PkgTemplates

t = Template(;
           user="arkinjo",
           authors=["Akira Kinjo"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )

t("EvoDevo2")
