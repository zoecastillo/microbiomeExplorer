/*taken from https://stackoverflow.com/questions/31703241/activate-tabpanel-from-another-tabpanel/31719425#31719425 */

shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value="' + name + '"]');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabledtab');
};

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value="' + name + '"]');
  tab.unbind('click.tab');
  tab.removeClass('disabledtab');
};

shinyjs.enableAll = function(){
  
  //$('.selectize-control .shinyjs-resettable').removeClass('disabled');
  //$('.selectize-dropdown').removeClass('disabled');
  $('div').removeClass('disabled');
};

shinyjs.collapse = function(boxid) {
  $('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
};