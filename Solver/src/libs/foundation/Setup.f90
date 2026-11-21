      Module mainKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: discretizationNodesKey     = "discretization nodes"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: inviscidDiscretizationKey  = "inviscid discretization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: splitFormKey               = "split form"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: meshFileNameKey            = "mesh file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: stlFileNameKey             = "stl file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartKey                 = "restart"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartFileNameKey         = "restart file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: solutionFileNameKey        = "solution file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveGradientsToSolutionKey = "save gradients with solution"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveLambVectorKey          = "save lamb vector"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveLambVectorDivergenceKey = "save lamb vector divergence"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveSensorToSolutionKey    = "save sensor with solution"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveLESToSolutionKey       = "save les with solution"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: saveSourceToSolutionKey    = "save source with solution"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfTimeStepsKey       = "number of time steps"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: outputIntervalKey          = "output interval"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: convergenceToleranceKey    = "convergence tolerance"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfBoundariesKey      = "number of boundaries"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(12) :: mainKeywords =  [ meshFileNameKey,           &
                                                                          inviscidDiscretizationKey,  &
                                                                          splitFormkey,               &
                                                                          discretizationNodesKey,     &
                                                                          saveGradientsToSolutionKey, &
                                                                          saveSensorToSolutionKey,    &
                                                                          restartKey,                 &
                                                                          restartFileNameKey,         &
                                                                          solutionFileNameKey,        &
                                                                          numberOfTimeStepsKey,       &
                                                                          outputIntervalKey,          &
                                                                          convergenceToleranceKey  ]
      END MODULE mainKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
