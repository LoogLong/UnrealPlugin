// Copyright LoogLong. All Rights Reserved.
#include "AnimGraphNode_SpringMagic.h"
#include "SpringMagicEditMode.h"
#include "Kismet2/CompilerResultsLog.h"

#define LOCTEXT_NAMESPACE "SpringMagic"

// ----------------------------------------------------------------------------
UAnimGraphNode_SpringMagic::UAnimGraphNode_SpringMagic(const FObjectInitializer& ObjectInitializer)
	: Super(ObjectInitializer)
{
}

FText UAnimGraphNode_SpringMagic::GetControllerDescription() const
{
	return LOCTEXT("SpringMagic", "Spring Magic");
}


// ----------------------------------------------------------------------------
FText UAnimGraphNode_SpringMagic::GetNodeTitle(ENodeTitleType::Type TitleType) const
{
	return LOCTEXT("NodeTitle", "Spring Magic");
}

void UAnimGraphNode_SpringMagic::PostEditChangeProperty(struct FPropertyChangedEvent& PropertyChangedEvent)
{
	Super::PostEditChangeProperty(PropertyChangedEvent);

}

FEditorModeID UAnimGraphNode_SpringMagic::GetEditorMode() const
{
	return FSpringMagicEditMode::ModeName;
}


void UAnimGraphNode_SpringMagic::Serialize(FArchive& Ar)
{
	Super::Serialize(Ar);
}
#undef LOCTEXT_NAMESPACE
