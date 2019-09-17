close session  
cd C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\  
set projection orthographic   
set bgColor white   
open TrajBox.bild  
scale 1.0 
windowsize 960 960 
~modeldisplay  #0 
 
open SSL_NoOH_Final_CG_1.bild  
turn x  models #1 180   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\image_1_.png  
turn y  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Sideimage_1_.png  
turn y  models #1 -90   
turn x  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Topimage_1_.png  
close #1 
 
open SSL_NoOH_Final_CG_2.bild  
turn x  models #1 180   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\image_2_.png  
turn y  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Sideimage_2_.png  
turn y  models #1 -90   
turn x  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Topimage_2_.png  
close #1 
 
open SSL_NoOH_Final_CG_3.bild  
turn x  models #1 180   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\image_3_.png  
turn y  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Sideimage_3_.png  
turn y  models #1 -90   
turn x  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Topimage_3_.png  
close #1 
 
open SSL_NoOH_Final_CG_4.bild  
turn x  models #1 180   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\image_4_.png  
turn y  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Sideimage_4_.png  
turn y  models #1 -90   
turn x  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Topimage_4_.png  
close #1 
 
open SSL_NoOH_Final_CG_5.bild  
turn x  models #1 180   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\image_5_.png  
turn y  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Sideimage_5_.png  
turn y  models #1 -90   
turn x  models #1 90   
copy file C:\Users\huang.2011\Dropbox\ShareFolderDNA\DOMM meeting\Draft\Released Code\SLL_wTraj\BILDs\Topimage_5_.png  
