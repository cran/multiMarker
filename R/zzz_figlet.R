.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive())
  { # colossal -- try also Nancyj
    packageStartupMessage(
      "
                                 M    M
      m    m u   u l     ttttt i MM  MM   aaa  rrrr  k   k eeeee rrrr
      mm  mm u   u l       t   i M MM M  a   a r   r k  k  e     r   r
      m mm m u   u l       t   i M    M  aaaaa rrrr  kkk   eee   rrrr
      m    m u   u l       t   i M    M  a   a r  r  k  k  e     r  r
      m    m uuuuu lllll   t   i M    M  a   a r   r k   k eeeee r   r

      ",
      "version ", version, "\n" )
  }
  else
  { packageStartupMessage("Package 'multiMarker' version ", version) }

  packageStartupMessage("Type 'citation(\"multiMarker\")' for citing this R package in publications.")
  invisible()
}
