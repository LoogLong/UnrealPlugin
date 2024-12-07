// Copyright LoogLong. All Rights Reserved.
#pragma once
#include "CoreMinimal.h"
#include "AnimNodeEditMode.h"
#include "AnimGraphNode_SpringMagic.h"
#include "AnimNode_SpringMagic.h"


class FEditorViewportClient;
class FPrimitiveDrawInterface;
class USkeletalMeshComponent;
struct FViewportClick;

class FSpringMagicEditMode : public FAnimNodeEditMode
{
public:
	FSpringMagicEditMode();

	/** IAnimNodeEditMode interface */
	virtual void EnterMode(class UAnimGraphNode_Base* InEditorNode, struct FAnimNode_Base* InRuntimeNode) override;
	virtual void ExitMode() override;

	/** FEdMode interface */
	virtual void Render(const FSceneView* View, FViewport* Viewport, FPrimitiveDrawInterface* PDI) override;
	virtual void DrawHUD(FEditorViewportClient* ViewportClient, FViewport* Viewport, const FSceneView* View,
	                     FCanvas* Canvas) override;
	static FEditorModeID ModeName;
private:
	
	/** Draw text func for DrawHUD */
	void DrawTextItem(const FText& Text, FCanvas* Canvas, float X, float& Y, float FontHeight);
	void Draw3DTextItem(const FText& Text, FCanvas* Canvas, const FSceneView* View, const FViewport* Viewport,
	                    FVector Location);

	/** Cache the typed nodes */
	struct FAnimNode_SpringMagic* RuntimeNode;
	UAnimGraphNode_SpringMagic* GraphNode;

};
