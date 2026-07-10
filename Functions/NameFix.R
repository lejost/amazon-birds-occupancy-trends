# --------------------------------------
# FUNCTION name_fix
# required packages: none
# description: takes species name from user input and makes sure it:
#          1. Corresponds to a species that appears in the SACC list
#          2. If it is a scientific name, translate to English
#          3. If it is one of the English names that are different between
#             SACC and BirdNet, convert to the BirdNet form
#          Needs to have the sacc names object in the workspace
# inputs: 'sps', a character string containing species names. 
# outputs: 'sps' output name
########################################

name_fix <- function(sps) {
  
  # buscar nomes de especies que aparecem diferente no SACC e no BirdNet
  odd_BNnames <- c("Yellow-margined Flycatcher","Gray-crowned Flycatcher",
                   "Guianan Trogon","Long-billed Gnatwren",
                   "Waved Woodpecker", "Long-tailed Woodcreeper",
                   "Plain Xenops", "Black-throated Trogon",
                   "Mealy Parrot", "Barn Owl", "Red-lored Parrot")
  
  odd_SACCnames <- c("Yellow-margined Flatbill","Gray-crowned Flatbill",
                     "Guianan Violaceous-Trogon","Trilling Gnatwren", 
                     "Variable Woodpecker", "Whistling Long-tailed Woodcreeper",
                     "Atlantic Plain-Xenops", "Amazonian Black-throated Trogon",
                     "Mealy Amazon", "American Barn Owl","Red-lored Amazon")

  # se for um nome ingles unico do BirdNet, converter temporariamente para SACC
  if(sps %in% odd_BNnames) { sps <- odd_SACCnames[which(odd_BNnames == sps)] }
  
  # se for um nome cientifico passar para ingles de acordo com SACC
  if(sps %in% sacc[,5]) { sps <- sacc[,6][which(sacc[,5]==sps)] }
  
  # se for um nome em ingles que nao existe no SACC retornar erro
  if(!sps %in% sacc[,6]) { stop("sps name not in SACC list") }
  
  # se for um nome ingles unico do SACC, converter novamente para nome do BirdNet
  # if(sps %in% odd_SACCnames) { sps <- odd_BNnames[which(odd_SACCnames == sps)] }
  
  return(sps)
  
}
