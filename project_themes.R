
# library(broman)
getMyProjectPalette <- function(){
  myProjectPalette = c(
    # match paper figures
    before = "#2E5473",
    'before opening'= "#2E5473",
    Preopening = "#2E5473",
    after = "#BB302F",
    'after opening'= "#BB302F",
    Postopening = "#BB302F",
    # Gloves (with extended categories)
    Glove1 = "deepskyblue",
    Glove2 = "deepskyblue3",
    'Glove from Glove Box' = "deepskyblue1",
    Glove = "deepskyblue1",
    'blank control' = "black", #white? consider: ("#fefe22") brocolors("crayons")["Laser Lemon"]
    # Staff_Surface >> muted blue
    Staff_Surface = "#a2a2d0", # brocolors("crayons")["Blue Bell"]
    'Personal Cell Phone' = "#a2a2d0", # brocolors("crayons")["Blue Bell"]
    'Hospital Pager' = "#dbd7d2", # brocolors("crayons")["Timberwolf"]
    'Shirt Hem' = "#c5d0e6", # brocolors("crayons")["Periwinkle"]
    Shoe = "#b0b7c6", # brocolors("crayons")["Cadet Blue"]
    # Station_Surface >> purples
    'Station_Surface' = "#f664af", # brocolors("crayons")["Magenta"]
    'Countertop' = "#dd4492", # brocolors("crayons")["Cerise"] 
    'Computer Mouse' = "#c364c5", # brocolors("crayons")["Fuchsia"]
    'Station Phone'= "#f664af", # brocolors("crayons")["Magenta"]
    'Chair Armrest' = "#de5d83", # brocolors("crayons")["Blush"]
    'Corridor Floor' = "purple", 
    # (Room_Surface | Station_Surface) > Water Faucet Handle >> green
    Room_Surface = "#71bc78", # brocolors("crayons")["Fern"]  
    'Hot Tap Water Faucet Handle' = "#71bc78", # brocolors("crayons")["Fern"]  
    'Cold Tap Water Faucet Handle'= "#87a96b", # brocolors("crayons")["Asparagus"]
    'Cold Tap Faucet Handle'="#a8e4a0" , # brocolors("crayons")["Granny Smith Apple"]
    # Patient_Skin | Staff_Skin >> organic?
    Staff_Skin = "#ffa089" , # brocolors("crayons")["Vivid Tangerine"]
    Patient_Skin = "#ffa089" , # brocolors("crayons")["Vivid Tangerine"]
    'Inguinal Fold'="#ffa089" , # brocolors("crayons")["Vivid Tangerine"]
    Hand="#fae7b5" , # brocolors("crayons")["Banana Mania"]
    Nose = "#ffcf48" , # brocolors("crayons")["Sunglow"]
    # Nose="#c5e384" , # brocolors("crayons")["Yellow Green"]
    Axilla="#f78fa7" , # brocolors("crayons")["Pink Sherbert"]
    # Room_Surface >> greenish/blueish
    Room_Surface = "#2b6cc4" , # brocolors("crayons")["Denim"]
    Floor = "#2b6cc4" , # brocolors("crayons")["Denim"]
    Bedrail = "#aaf0d1", # brocolors("crayons")["Magic Mint"]
    Water = "dodgerblue"
  )
  return(myProjectPalette)
}

# to use this with ggplot2, add a scale_*_manual layer:
# plot1 + 
#    scale_colour_manual(values = myProjectPalette) 