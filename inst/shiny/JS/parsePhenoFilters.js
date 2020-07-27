
shinyjs.parsePhenoFilters = 
  function(){
  var phenorows = document.getElementsByClassName('phenoremoverow');
  var phenoids = [];
  l = phenorows.length;
    for (i = 0; i < l; i++) {
      phenoids.push(phenorows[i].id);
    }
  Shiny.setInputValue('loadnfilter-parsedFilters',phenoids);
};