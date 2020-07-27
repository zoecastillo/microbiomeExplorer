
shinyjs.getTraceName = 
  function(params){
  var plotdata = params[0];
  var curvenum = params[1];
  if(plotdata !== null){
    Shiny.setInputValue('intraAnalysis-abundancePlot-clickedFeature',plotdata['x']['data'][curvenum]['name']);
  }
};