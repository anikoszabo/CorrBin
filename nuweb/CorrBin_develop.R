library(devtools)
source('nuweb/Nuweb.R')

cb <- as.package(".")
nuweb(cb)


# create native routine registration file
tools::package_native_routine_registration_skeleton("z://RForge/corrbin/pkg/CorrBin",
                                                    con="z://RForge/corrbin/pkg/CorrBin/src/init.c")

document()
run_examples()  # or dev_example("ran.CMData")
load_all()

build()
check(check_dir = "c:/TEMP", cran = TRUE)
