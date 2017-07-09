

require('shiny')
if(!require('covarShiny')){
  install.packages("./covarShiny/",repos=NULL, type='source' )
  library('covarShiny')
}
library('covarShiny')

#add ui elementrequire()s as Input() and Output() functions as arguments into fluidpage
ui=fluidPage(
  titlePanel('Enter FASTA sequences for covariation analysis'),
  
  sidebarLayout(
    #inputs  
    sidebarPanel(
      fileInput(inputId='f',label='Input a file containing FASTA sequences.'), 
      numericInput(inputId='pos1', label='position 1', value=NA, min=1, max= NA,step= 100)    ,
      numericInput(inputId='pos2', label='position 2',value=NA, min=1, max= NA,step= 100))    ,
    #end sidebar
    mainPanel(
      tabsetPanel(
        tabPanel("table", tableOutput('matrix')),
        tabPanel("map", imageOutput('image'))
      )
    )
  ))





server=function(input,output){
  {
    
    

      output$matrix = renderTable({inFile=input$f
                                        covImageui(covarRNAui(inFile$datapath,input$pos1,input$pos2),filename='heatmap.jpeg')})
      output$image=renderImage({
        
        input$f
        input$pos1
        input$pos2
       list(src='./covarImages/heatmap.jpeg',  contentType = 'image/jpeg' ,width=500, height=500)
      
    })  
    
    
    
  }
  
  
}

shinyApp(server=server,ui=ui)

