try import DrWatson
catch
    import Pkg
    Pkg.add("DrWatson")
end
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Rath2017")

## ------------------------------------------------------------------------
# Download raw

## ------------------------------------------------------------------------
# Install unregistered packages
using Pkg
pkg"rm Chemostat"
pkg"rm UtilsJL"
pkg"add https://github.com/josePereiro/UtilsJL.git#v0.2.3"
pkg"add https://github.com/josePereiro/Chemostat#v0.7.0"
pkg"instantiate"
pkg"build"