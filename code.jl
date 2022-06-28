using PkgTemplates

tpl = Template(;
    user = "ManuelStapper",
    authors = "Manuel Stapper",
    dir="~/code",
    plugins=[
        Git(; manifest=true, ssh=true),
        Codecov(),
    ],
)

cd("C:\\Users\\stapperm\\.julia\\dev\\RandomVariables")
tpl("RandomVariables.jl")
