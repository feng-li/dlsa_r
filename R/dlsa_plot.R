#' Plot for the results
#'
#' @param plot_fun user-submitted plotting function
#' @param roc logical. Should the roc curve be included in the final result?
#' @param Y a factor, numeric or character vector of responses, typically encoded with 0 (controls) and 1 (cases)
#' @param p a numeric or ordered vector of the same length than response, containing the predicted probability of each observation. 
#' @param ... optional arguments to plot_fun
#' @return A list of the plots
#' @importFrom pROC roc auc ggroc
#' @import ggplot2
#' @export
#' @examples
#' # Not Run
dlsa.plot <- function(plot_fun, roc = 0, Y = NULL, p = NULL, ...)
{
  if(roc == 1){
    
    if(is.null(Y) & is.null(p)) stop('please input the arguments for roc plot')
    
    dlsa_roc = roc(Y, p)
    dlsa_roc_plot = ggroc(dlsa_roc,size = 1)+theme_minimal() + ggtitle("DLSA ROC curve") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
    dlsa_auc = auc(dlsa_roc)
    
    plot = plot_fun(...)
    
    return(list(roc = dlsa_roc_plot, auc = dlsa_auc, plot = plot))
  }else{
    
    plot = plot_fun(...)
    return(list(plot = plot))
    
  }
}
  


