# Part of the mpt2irt package for estimating model parameters
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
    mpt2irtLib <- dirname(system.file(package = "mpt2irt"))
    pkgdesc <- suppressWarnings(utils::packageDescription("mpt2irt", lib.loc = mpt2irtLib))
    if (length(pkgdesc) > 1) {
        # builddate <- gsub(';.*$', '', pkgdesc$Built)
        packageStartupMessage(paste("mpt2irt (Version ", pkgdesc$Version,
                                    # ", packaged: ", builddate,
                                    ")", sep = ""))
    }
    # packageStartupMessage("- Do not expect the default priors to remain the same in future mpt2irt versions.")
    # packageStartupMessage("Thus, R scripts should specify priors explicitly, even if they are just the defaults.")
    # packageStartupMessage("- For execution on a local, multicore CPU with excess RAM we recommend calling")
    # packageStartupMessage("options(mc.cores = parallel::detectCores())")
}
