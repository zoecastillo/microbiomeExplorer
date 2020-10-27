

shinyjs.resetAxes = function(){ 
  var reset_button = $('#interAnalysis-betaDiv-betaDiv .modebar-btn');
  if (typeof reset_button !== 'undefined'){
    if(reset_button.length >= 5){
      reset_button[5].click();
    }
  }
};