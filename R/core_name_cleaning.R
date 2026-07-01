core_name_cleaning <- function(name, empty_place_holder="")
{

  name <- toupper(name)
  # replace < <= > >= and ==
  name <- gsub("!=","_diseq_",name)
  name <- gsub(">=","_gte_",name)
  name <- gsub("> =","_gte_",name)
  name <- gsub("=>","_gte_",name)
  name <- gsub("= >","_gte_",name)
  name <- gsub("<=","_lte_",name)
  name <- gsub("< =","_lte_",name)
  name <- gsub("=<","_lte_",name)
  name <- gsub("= <","_lte_",name)
  name <- gsub("==","_eq_",name)
  name <- gsub("=", "_eq_", name)

  name <- gsub("<","_lt_",name)
  name <- gsub(">","_gt_",name)

  name <- toupper(gsub("[.]", "_", name))
  name <- toupper(gsub(" ", "_", name))
  name <- toupper(gsub("-", "_", name))
  name <- toupper(gsub("'", "_", name))
  name <- toupper(gsub("__", "_", name))
  name <- toupper(gsub("'", "_", name))


  # replace any special character with underscore
  name <- gsub("[^a-zA-Z0-9]", "_", name)
  name <- toupper(gsub("__", "_", name))
  # remove any underscore at the end
  name <- sub("_+$", "", name)

  # remove any underscore at the beginning
  name <- sub("^_+", "", name)

  # remove any dot at the end
  name <- sub("[.]+$", "", name)

  # remove any dot at the beginning
  name <- sub("^[.]+", "", name)
  name <- toupper(gsub("__", "_", name))
  name <- toupper(gsub("__", "_", name))
  name <- toupper(gsub("__", "_", name))

  if (any(name=="") & (empty_place_holder!=""))
      name[name==""] <- empty_place_holder

  return(name)
}
