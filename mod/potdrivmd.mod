	  V(  d   k820309              15.0        qQf                                                                                                           
       src/potdriv.f POTDRIVMD                                                     
       TASKNUM                                                          #         @                                                        #         @                                                   
   #ZZ    #BETANH3    #ALPHAS    #R1    #BETAS 	   #GAMMAS 
   #ROW1    #ROW2    #THETA    #VX              
  @                                   
                
  @                                   
                
  @                                   
                
  @                                   
                
  @                              	     
                
  @                              
     
                
  @                                   
                
  @                                   
                
  @                                   
                D @                                   
       #         @                                                  
   #CART_HNH3_TO_RADAU%ABS    #CART_HNH3_TO_RADAU%SIN    #CART_HNH3_TO_RADAU%ACOS    #CART_HNH3_TO_RADAU%SQRT    #CART_HNH3_TO_RADAU%SUM    #CART0    #RS    #BETA_NH3    #ALPHA_S    #R1    #BETA_S    #GAMMA_S    #ROW1    #ROW2    #THETA                                                   ABS                                                SIN                                                ACOS                                                SQRT                                                SUM           
                                                    
    p          p          p            p          p                                    D                                     
                 D @                                   
                 D                                     
                 D                                     
                 D @                                   
                 D                                     
                 D @                                   
                 D @                                   
                 D @                                   
       #         @                                                       #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_BOTTOM     #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_IN_PLACE "   #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ARGV_NULL $   #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ARGVS_NULL &   #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ERRCODES_IGNORE (   #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_STATUS_IGNORE *   #POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_STATUSES_IGNORE ,                                                                        #NH4_POT_ENG_SUF%MPI_FORTRAN_BOTTOM%MPI_BOTTOM !             �            �                  !                                                              "                          #NH4_POT_ENG_SUF%MPI_FORTRAN_IN_PLACE%MPI_IN_PLACE #             �            �                  #                                                              $                          #NH4_POT_ENG_SUF%MPI_FORTRAN_ARGV_NULL%MPI_ARGV_NULL %   -          �            �                  %                                 p          p            p                                                                            &                          #NH4_POT_ENG_SUF%MPI_FORTRAN_ARGVS_NULL%MPI_ARGVS_NULL '             �            �                  '             
                                                 (                          #NH4_POT_ENG_SUF%MPI_FORTRAN_ERRCODES_IGNORE%MPI_ERRCODES_IGNORE )             �            �                  )                                 p          p            p                                                                            *                          #NH4_POT_ENG_SUF%MPI_FORTRAN_STATUS_IGNORE%MPI_STATUS_IGNORE +             �            �                  +                                 p          p            p                                                                            ,                          #NH4_POT_ENG_SUF%MPI_FORTRAN_STATUSES_IGNORE%MPI_STATUSES_IGNORE -             �            �                  -             
       #         @                                  .                
   #RADAU_HNH3_TO_CART%SUM /   #RS 0   #BETA_NH3 1   #ALPHA_S 2   #R1 3   #BETA_S 4   #GAMMA_S 5   #ROW1 6   #ROW2 7   #THETA 8   #CART0 9                                             /     SUM           
                                 0     
                
@ @                              1     
                
@ @                              2     
                
                                 3     
                
@ @                              4     
                
@ @                              5     
                
  @                              6     
                
  @                              7     
                
  @                              8     
                D @                              9                   
     p          p          p            p          p                          #         @                                  :                    #CARTIN ;   #CARTOUT <             
  @                              ;                   
 %   p          p          p            p          p                                    D @                              <                   
 &    p          p          p            p          p                          #         @                                  =                   #RADAU_CART%COS >   #RADAU_CART%SIN ?   #RADAU_CART%SQRT @   #RADAU_CART%SUM A   #AMASS B   #ROW1 C   #ROW2 D   #THETA E   #CART0 F                                             >     COS                                           ?     SIN                                           @     SQRT                                           A     SUM           
                                 B                   
    p          p            p                                    
                                 C     
                
                                 D     
                
  @                              E     
                D                                F     	              
     p          p          p            p          p                          #         @                                  G                   #ROTEULAR%DCOS H   #ROTEULAR%DSIN I   #NATOM J   #CARTES K   #A L   #B M   #C N                                             H     DCOS                                           I     DSIN           
                                  J                    
D                                K                    
 #    p          p          5 � p        r J       p          5 � p        r J                                @                              L     
                  @                              M     
                  @                              N     
       #         @                                  O                   #CART_RADAU%ACOS P   #CART_RADAU%SQRT Q   #CART_RADAU%SUM R   #AMASS S   #CART0 T   #ROW1 U   #ROW2 V   #THETA W   #AP X                                             P     ACOS                                           Q     SQRT                                           R     SUM           
                                 S                   
    p          p            p                                    
                                 T     	              
    p          p          p            p          p                                    D                                U     
                 D                                V     
                 D @                              W     
                 D                                X                   
     p          p            p                          #         @                                  Y                    #NDIM Z   #XIN [   #XOUT \             
                                  Z                    
                                 [                    
 '   p          5 � p        r Z       5 � p        r Z                              D                                \                    
 (    p          5 � p        r Z       5 � p        r Z                     #         @                                   ]                   #CART_BOND%DACOS ^   #CART_BOND%ACOS _   #CART_BOND%SQRT `   #CART_BOND%SUM a   #CART b   #ROUT c                                             ^     DACOS                                           _     ACOS                                           `     SQRT                                           a     SUM           
                                 b                   
 )   p          p          p            p          p                                    D                                c                   
 *    p          p            p                             �          fn#fn    �   H   J  DECDATAMD "     @       TASKNUM+DECDATAMD    H  H       POT_ZERO    �  �       POTDRIV    ?  @   a   POTDRIV%ZZ       @   a   POTDRIV%BETANH3    �  @   a   POTDRIV%ALPHAS    �  @   a   POTDRIV%R1    ?  @   a   POTDRIV%BETAS      @   a   POTDRIV%GAMMAS    �  @   a   POTDRIV%ROW1    �  @   a   POTDRIV%ROW2    ?  @   a   POTDRIV%THETA      @   a   POTDRIV%VX #   �  D      CART_HNH3_TO_RADAU '     <      CART_HNH3_TO_RADAU%ABS '   ?  <      CART_HNH3_TO_RADAU%SIN (   {  =      CART_HNH3_TO_RADAU%ACOS (   �  =      CART_HNH3_TO_RADAU%SQRT '   �  <      CART_HNH3_TO_RADAU%SUM )   1  �   a   CART_HNH3_TO_RADAU%CART0 &   �  @   a   CART_HNH3_TO_RADAU%RS ,   %  @   a   CART_HNH3_TO_RADAU%BETA_NH3 +   e  @   a   CART_HNH3_TO_RADAU%ALPHA_S &   �  @   a   CART_HNH3_TO_RADAU%R1 *   �  @   a   CART_HNH3_TO_RADAU%BETA_S +   %	  @   a   CART_HNH3_TO_RADAU%GAMMA_S (   e	  @   a   CART_HNH3_TO_RADAU%ROW1 (   �	  @   a   CART_HNH3_TO_RADAU%ROW2 )   �	  @   a   CART_HNH3_TO_RADAU%THETA     %
  �      NH4_POT_ENG_SUF =   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_BOTTOM >   p  H      NH4_POT_ENG_SUF%MPI_FORTRAN_BOTTOM%MPI_BOTTOM ?   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_IN_PLACE B   ?  H      NH4_POT_ENG_SUF%MPI_FORTRAN_IN_PLACE%MPI_IN_PLACE @   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ARGV_NULL D     �      NH4_POT_ENG_SUF%MPI_FORTRAN_ARGV_NULL%MPI_ARGV_NULL A   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ARGVS_NULL F   ?  H      NH4_POT_ENG_SUF%MPI_FORTRAN_ARGVS_NULL%MPI_ARGVS_NULL F   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_ERRCODES_IGNORE P     �      NH4_POT_ENG_SUF%MPI_FORTRAN_ERRCODES_IGNORE%MPI_ERRCODES_IGNORE D   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_STATUS_IGNORE L   Q  �      NH4_POT_ENG_SUF%MPI_FORTRAN_STATUS_IGNORE%MPI_STATUS_IGNORE F   �  �   �   POTDRIVMD!NH4_POT_ENG_SUF!MPI_FORTRAN_STATUSES_IGNORE P   �  H      NH4_POT_ENG_SUF%MPI_FORTRAN_STATUSES_IGNORE%MPI_STATUSES_IGNORE #   �  �       RADAU_HNH3_TO_CART '   �  <      RADAU_HNH3_TO_CART%SUM &   �  @   a   RADAU_HNH3_TO_CART%RS ,      @   a   RADAU_HNH3_TO_CART%BETA_NH3 +   `  @   a   RADAU_HNH3_TO_CART%ALPHA_S &   �  @   a   RADAU_HNH3_TO_CART%R1 *   �  @   a   RADAU_HNH3_TO_CART%BETA_S +      @   a   RADAU_HNH3_TO_CART%GAMMA_S (   `  @   a   RADAU_HNH3_TO_CART%ROW1 (   �  @   a   RADAU_HNH3_TO_CART%ROW2 )   �  @   a   RADAU_HNH3_TO_CART%THETA )      �   a   RADAU_HNH3_TO_CART%CART0    �  a       REARRANGE_ATOM &   5  �   a   REARRANGE_ATOM%CARTIN '   �  �   a   REARRANGE_ATOM%CARTOUT    �  �       RADAU_CART    k  <      RADAU_CART%COS    �  <      RADAU_CART%SIN     �  =      RADAU_CART%SQRT       <      RADAU_CART%SUM !   \  �   a   RADAU_CART%AMASS     �  @   a   RADAU_CART%ROW1     0  @   a   RADAU_CART%ROW2 !   p  @   a   RADAU_CART%THETA !   �  �   a   RADAU_CART%CART0    d  �       ROTEULAR    �  =      ROTEULAR%DCOS    ;  =      ROTEULAR%DSIN    x  @   a   ROTEULAR%NATOM     �  �   a   ROTEULAR%CARTES    �  @   a   ROTEULAR%A    �  @   a   ROTEULAR%B      @   a   ROTEULAR%C    L  �       CART_RADAU        =      CART_RADAU%ACOS     L   =      CART_RADAU%SQRT    �   <      CART_RADAU%SUM !   �   �   a   CART_RADAU%AMASS !   Y!  �   a   CART_RADAU%CART0     "  @   a   CART_RADAU%ROW1     M"  @   a   CART_RADAU%ROW2 !   �"  @   a   CART_RADAU%THETA    �"  �   a   CART_RADAU%AP    a#  e       DCOPY1D    �#  @   a   DCOPY1D%NDIM    $  �   a   DCOPY1D%XIN    �$  �   a   DCOPY1D%XOUT    n%  �       CART_BOND     &  >      CART_BOND%DACOS    X&  =      CART_BOND%ACOS    �&  =      CART_BOND%SQRT    �&  <      CART_BOND%SUM    '  �   a   CART_BOND%CART    �'  �   a   CART_BOND%ROUT 