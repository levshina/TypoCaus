#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#Author: Natalia Levshina
#Version: 28.06.2018

database <- read.table("typocaus_main.txt", header = T, quote = "", sep = "\t", comment = "", encoding = "UTF-8")

library(shiny)
library(networkD3)
library(ca)
library(smacof)
library(cluster)

semantics <- database[, c(20:34)]
semantics <- as.matrix(semantics)
semantics[semantics == "No"] <- NA



# Define UI 
ui <- fluidPage(
   
   # Application title
   titlePanel("Database of causative constructions in languages of the world"),
   tabsetPanel(type = "tabs",

tabPanel("Search by language", 
   sidebarLayout(
     
     # Sidebar panel for inputs ----
     sidebarPanel(
       
       # Input: Select a dataset ----
       textInput("languageName", "Enter the name of a language:" ),
       #hr(),
       #textInput("ISO3", "Enter the 3-digit code of a language:" ),


       # Input: actionButton() to defer the rendering of output ----
       # until the user explicitly clicks the button (rather than
       # doing it immediately when inputs change). This is useful if
       # the computations required to render output are inordinately
       # time-consuming.
       actionButton("updateLanguage", "Submit")
     ),
     
     # Main panel for displaying outputs ----
     mainPanel(
       
       # Output: Header + summary of distribution ----
       verbatimTextOutput("metainformation"),
       hr(),
       h4(textOutput("cxnnumber")),
       htmlOutput("constructions")
            )
   )), 


tabPanel("Semantic maps", 
         sidebarLayout(
           
           # Sidebar panel for inputs ----
           sidebarPanel(
             
             radioButtons("maptype", label = h3("Type of semantic map"),
                          choices = list("Network of functions" = 'network', "Distances between functions" = 'mds', "Distances between constructions and between functions" = 'ca'), 
                          selected = 'network'),
             hr(),
             checkboxGroupInput("show_vars", "Semantic functions to include:",
                                colnames(semantics), selected = colnames(semantics)),
             actionButton("updateMap", "Submit")
             
             ),
           
           # Main panel for displaying outputs ----
           mainPanel(
             conditionalPanel(
               condition = "input.maptype == 'network'",
               forceNetworkOutput("semnetwork")
               ),
             conditionalPanel(
               condition = "input.maptype == 'mds'",
               plotOutput("mds_map")
             ),
             conditionalPanel(
               condition = "input.maptype == 'ca'",
               plotOutput("ca_map")
             )
             
             
         #forceNetworkOutput("semnetwork"))
)
)
)))



server <- function(input, output) {

  languageInput <- eventReactive(input$updateLanguage,
                                 {lang <- input$languageName},
                                 ignoreNULL = TRUE)

  output$metainformation <- renderPrint({
    lang <- languageInput()
    print (lang)
    print(unique(database[database$Languoid == lang, 2:6]), row.names  = F)
  })  
  

  output$cxnnumber <- renderText({
    lang <- languageInput()
    langcaus <- database[database$Languoid == lang,]
    total <- as.numeric(nrow(langcaus))
    print(paste(as.character(total), " causative construction(s) found!"), quote = FALSE, row.names = F)
  })
  
  output$constructions <- renderUI({
    lang <- languageInput()
    langcaus <- database[database$Languoid == lang,]
    output <- "" 
    for (i in 1:nrow(langcaus)){
    output_id <- h4(paste("Construction ", i, ": " , as.character(langcaus[i, 12])))
    output_form <- paste("<b>Form: </b>", langcaus[i, 7])
    output_meaning <- paste("<b>Meaning: </b>", langcaus[i, 8])
    output_example <- ifelse(is.na(langcaus[i, 10]), as.character(langcaus[i, 9]), paste(langcaus[i, 9], langcaus[i, 11], langcaus[i, 10], sep = "<br/>"))
    output_example1 <- paste("<b>Example: </b><br/>", output_example)
    newoutput <- paste(output_id, output_form, output_meaning, output_example1, sep = "<br/>")
    output <- paste(output, newoutput, sep = "<br/><br/>") 
    } 
    HTML(output)
  })

  
  
  
  mapInput <- eventReactive(input$updateMap,
                                 {type <- input$maptype})
  
  varInput <- eventReactive(input$updateMap,
                    {vars <- input$show_vars})
  
  output$semnetwork <- renderForceNetwork({
    cooccurrences <- matrix(0, ncol = 3)
    semantics1 <- semantics[, varInput()]
      for (j in 1:(ncol(semantics1)-1)){
        for (k in (j+1):ncol(semantics1)){
          if (j != k){  
            mytable <- table(semantics1[, j], semantics1[,k], useNA = "no")
            for (r in 1:nrow(mytable)){
              for (col in 1:ncol(mytable)){
                if (mytable[r, col] > 0){
                  cooccurrences <- rbind(cooccurrences, c(j, k, mytable[r, col]))  
                }  
              }  
            }
          }
        }
      }
    
    NetworkLinks <- as.data.frame(cooccurrences[-1, ])
    colnames(NetworkLinks) <- c("source", "target", "value")
    NetworkLinks$source = NetworkLinks$source -1
    NetworkLinks$target = NetworkLinks$target -1
    NetworkNodes <- data.frame(name = colnames(semantics1), group = rep(1, ncol(semantics1)), size = rep(30, ncol(semantics1)))  
  forceNetwork(Links = NetworkLinks, Nodes = NetworkNodes, Source = "source",
                 Target = "target", Value = "value", NodeID = "name",
                 Group = "group", opacity = 0.9, fontSize = 30, Nodesize = "size", 
               linkDistance = 150, linkWidth = JS("function(d) { return Math.sqrt(d.value)*2.5; }"), linkColour = "green", height = 20000, width = 10000, opacityNoHover = 0.8, bounded = F)
  })     
  
  output$mds_map <- renderPlot({
  semantics1 <- semantics[, varInput()] 
  semantics1[is.na(semantics1)] <- "No"
  semdist <- daisy(as.data.frame(t(semantics1)))
  semmds <- mds(semdist)
  plot(semmds$conf, type = "n", xlim = c(min(semmds$conf[, 1]) - 0.3, max(semmds$conf[, 1]) + 0.3), ylim = c(min(semmds$conf[, 2]) - 0.1, max(semmds$conf[, 2]) + 0.1))
  text(semmds$conf, labels = rownames(semmds$conf), cex = 1.2)
  })
  
  output$ca_map <- renderPlot({
  semantics1 <- semantics[, varInput()]
  contingency <- matrix(0, nrow = length(levels(database$Type_cat)), ncol = ncol(semantics1))
  for (j in 1:ncol(semantics1)){
  mytab <- table(database$Type_cat, semantics1[, j])  
  contingency[, j] <- mytab
  } 
  rownames(contingency) <- levels(database$Type_cat)
  colnames(contingency) <- colnames(semantics1)
  camap <- ca(contingency)
  par(cex = 1.5)
  plot(camap)  
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


