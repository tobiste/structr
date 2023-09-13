# shmax_from_tensor <- function(x, tensor = c("stress", "FMS")){
#   tensor <- match.arg(tensor)
#   
#   if(tensor == "stress"){
#   # Calculating the direction of maximum horizontal stress using the complete 
#   # stress tensor
#   dividend = 2 *(S1 * s1n * s1e + S2 * s2n * s2e + S3 * s3n * s3e)
#   divisor = S1 * (s1n^2 - s1e^2) + S2 * (s2n^2 - s2e^2) + S3 * (s3n^2 - s3e^2)
#   
#   tan_2a = dividend/ divisor
#   a = atand(tan_2a) / 2
#   
#   } else {
#     # Calculating the direction of maximum horizontal stress from focal 
#     # mechanism stress inversion results
#     
#     dividend = 2 *(S1 * s1n * s1e + (1 - R) * s2n * s2e)
#     divisor = (s1n^2 - s1e^2) + (1 - R) * (s2n^2 - s2e^2)
#     if(dividend == 0){
#      if(abs(s1n) = abs(s1e) & R==1){
#        if(){
#          SH = NA
#        } else {
#          ah = a1
#        }
#      } else if(abs(s1n) = abs(s1e) & abs(s2n) == abs(s2e)){
#        if(R == 0){
#          SH = NA
#        } else {
#          ah = a1
#        }
#      } else if(divisor == 0){
#        if(R == 0){
#          SH = NA
#        } else {
#          ah = a1
#        }
#      }
#     }
#   }
# }