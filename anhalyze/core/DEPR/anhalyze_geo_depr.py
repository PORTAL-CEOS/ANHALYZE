import deprecation


# TODO can make these sections into its own class
#  @deprecation.deprecated()
def getIndex_sec(sectName):
    """Return the number associated to each section.
    One can check the respective numbers and sections on the Lab Guide.

    Inputs
    sectName = Section name; type: String

    Options:
    Bering Strait
    Lancaster Sound
    Jones Sound
    Nares Strait
    Davis Strait
    Fram Strait

"""

# Todo, add the tilted sections. So far, the section functions only work with
# zonal or meridional oriented sections.

    if sectName == 'Bering Strait':
        sect = 3180
    elif sectName=='Lancaster Sound':
        sect = 4180
    elif sectName=='Jones Sound':
        sect = 4180
    elif sectName=='Nares Strait':
        sect = 4270
    elif sectName=='Davis Strait':
        sect = 9360
    elif sectName=='Fram Strait':
        sect = 5360
    else:
        sect = None
        print("Section name not found.\nCall help/doc to check the sections available in this function.")

    return sect
