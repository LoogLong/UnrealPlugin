// Copyright LoogLong. All Rights Reserved.
#include "SpringMagicEditMode.h"
#include "CanvasItem.h"
#include "CanvasTypes.h"
#include "EditorViewportClient.h"
#include "IPersonaPreviewScene.h"
#include "SceneManagement.h"
#include "Animation/DebugSkelMeshComponent.h"
#include "Materials/MaterialInstanceDynamic.h"

#define LOCTEXT_NAMESPACE "SpringMagicEditMode"


FEditorModeID FSpringMagicEditMode::ModeName = "AnimGraph.SkeletalControl.SpringMagic";
FSpringMagicEditMode::FSpringMagicEditMode()
	: RuntimeNode(nullptr)
	  , GraphNode(nullptr)
{
}

void FSpringMagicEditMode::EnterMode(UAnimGraphNode_Base* InEditorNode, FAnimNode_Base* InRuntimeNode)
{
	RuntimeNode = static_cast<FAnimNode_SpringMagic*>(InRuntimeNode);
	GraphNode = CastChecked<UAnimGraphNode_SpringMagic>(InEditorNode);

	FAnimNodeEditMode::EnterMode(InEditorNode, InRuntimeNode);
}

void FSpringMagicEditMode::ExitMode()
{
	GraphNode = nullptr;
	RuntimeNode = nullptr;

	FAnimNodeEditMode::ExitMode();
}

void FSpringMagicEditMode::Render(const FSceneView* View, FViewport* Viewport, FPrimitiveDrawInterface* PDI)
{
	const USkeletalMeshComponent* SkelMeshComp = GetAnimPreviewScene().GetPreviewMeshComponent();

	if (RuntimeNode != nullptr)
	{
		const auto& Particles = RuntimeNode->SimulationParticles;
		for (auto& JointChain : RuntimeNode->SimulationJoints)
		{
			for (FSpringMagicJoint& Chain : JointChain)
			{
				const auto& Grandparent = Particles[Chain.GrandParent];
				const auto& Parent = Particles[Chain.Parent];
				const auto& Child = Particles[Chain.Child];
				if (Grandparent.bRootBone)
				{
					PDI->DrawLine(Grandparent.SimulationTransform.GetTranslation(), Parent.SimulationTransform.GetTranslation(), FLinearColor::Black, SDPG_Foreground);
					DrawWireSphere(PDI, Grandparent.SimulationTransform, FLinearColor::Blue, 1.0f, 16, SDPG_Foreground);
					DrawWireSphere(PDI, Parent.SimulationTransform, FLinearColor::Black, 1.0f, 16, SDPG_Foreground);
				}
				PDI->DrawLine(Parent.SimulationTransform.GetTranslation(), Child.SimulationTransform.GetTranslation(), FLinearColor::Black, SDPG_Foreground);

				DrawWireSphere(PDI, Child.SimulationTransform, FLinearColor::Black, 1.0f, 16, SDPG_Foreground);
				DrawWireSphere(PDI, Chain.NewChildPosition, FLinearColor::Yellow, 1.2f, 16, SDPG_Foreground);
				DrawWireSphere(PDI, Chain.CurrentChildPosition, FLinearColor::White, 1.2f, 16, SDPG_Foreground);
			}
		}

		FMaterialRenderProxy* Material = GEngine->ConstraintLimitMaterialPrismatic->GetRenderProxy();
		for (auto& Sphere : RuntimeNode->SphereColliders)
		{
			DrawSphere(PDI, Sphere.CenterSimSpace, FRotator::ZeroRotator, FVector(Sphere.Radius), 16, 4, Material, SDPG_World);
		}
		for (auto& Capsule : RuntimeNode->CapsuleColliders)
		{
			DrawCylinder(PDI, Capsule.CenterASimSpace, Capsule.CenterBSimSpace, Capsule.Radius, 25, Material, SDPG_World);
			DrawSphere(PDI, Capsule.CenterASimSpace, FRotator::ZeroRotator, FVector(Capsule.Radius), 24, 6, Material, SDPG_World);
			DrawSphere(PDI, Capsule.CenterBSimSpace, FRotator::ZeroRotator, FVector(Capsule.Radius), 24, 6, Material, SDPG_World);
		}
	}
	FAnimNodeEditMode::Render(View, Viewport, PDI);
}

void FSpringMagicEditMode::DrawHUD(FEditorViewportClient* ViewportClient, FViewport* Viewport, const FSceneView* View,
                                     FCanvas* Canvas)
{
	// float FontWidth, FontHeight;
	// GEngine->GetSmallFont()->GetCharSize(TEXT('L'), FontWidth, FontHeight);
	// constexpr float XOffset = 5.0f;
	// float DrawPositionY = Viewport->GetSizeXY().Y / Canvas->GetDPIScale() - (3 + FontHeight) - 100 / Canvas->
	// 	GetDPIScale();
	//
	// if (!FAnimWeight::IsRelevant(RuntimeNode->GetAlpha()) || !RuntimeNode->IsRecentlyEvaluated())
	// {
	// 	DrawTextItem(
	// 		LOCTEXT("", "This node does not evaluate recently."), Canvas, XOffset, DrawPositionY,
	// 		FontHeight);
	// 	FAnimNodeEditMode::DrawHUD(ViewportClient, Viewport, View, Canvas);
	// 	return;
	// }
	//
	// DrawTextItem(LOCTEXT("", "Q : Cycle Transform Coordinate System"), Canvas, XOffset, DrawPositionY, FontHeight);
	// DrawTextItem(
	// 	LOCTEXT("", "Space : Cycle Between Translate, Rotate and Scale"), Canvas, XOffset, DrawPositionY, FontHeight);
	// DrawTextItem(LOCTEXT("", "R : Scale Mode"), Canvas, XOffset, DrawPositionY, FontHeight);
	// DrawTextItem(LOCTEXT("", "E : Rotate Mode"), Canvas, XOffset, DrawPositionY, FontHeight);
	// DrawTextItem(LOCTEXT("", "W : Translate Mode"), Canvas, XOffset, DrawPositionY, FontHeight);
	// DrawTextItem(LOCTEXT("", "------------------"), Canvas, XOffset, DrawPositionY, FontHeight);
	//
	//
	// FString CollisionDebugInfo = FString(TEXT("Select Collision : "));
	// switch (SelectCollisionType)
	// {
	// case ECollisionLimitType::Spherical:
	// 	CollisionDebugInfo.Append(FString(TEXT("Spherical")));
	// 	break;
	// case ECollisionLimitType::Capsule:
	// 	CollisionDebugInfo.Append(FString(TEXT("Capsule")));
	// 	break;
	// case ECollisionLimitType::Box:
	// 	CollisionDebugInfo.Append(FString(TEXT("Box")));
	// 	break;
	// case ECollisionLimitType::Planar:
	// 	CollisionDebugInfo.Append(FString(TEXT("Planar")));
	// 	break;
	// default:
	// 	CollisionDebugInfo.Append(FString(TEXT("None")));
	// 	break;
	// }
	// if (SelectCollisionIndex >= 0)
	// {
	// 	CollisionDebugInfo.Append(FString(TEXT("[")));
	// 	CollisionDebugInfo.Append(FString::FromInt(SelectCollisionIndex));
	// 	CollisionDebugInfo.Append(FString(TEXT("]")));
	// }
	// DrawTextItem(FText::FromString(CollisionDebugInfo), Canvas, XOffset, DrawPositionY, FontHeight);
	//
	// const UDebugSkelMeshComponent* PreviewMeshComponent = GetAnimPreviewScene().GetPreviewMeshComponent();
	// if (GraphNode->bEnableDebugBoneLengthRate)
	// {
	// 	if (PreviewMeshComponent != nullptr && PreviewMeshComponent->MeshObject != nullptr)
	// 	{
	// 		for (auto& Bone : RuntimeNode->ModifyBones)
	// 		{
	// 			// Refer to FAnimationViewportClient::ShowBoneNames
	// 			const FVector BonePos = PreviewMeshComponent->GetComponentTransform().TransformPosition(Bone.Location);
	// 			Draw3DTextItem(FText::AsNumber(Bone.LengthRateFromRoot), Canvas, View,
	// 			               Viewport, BonePos);
	// 		}
	// 	}
	// }

	FAnimNodeEditMode::DrawHUD(ViewportClient, Viewport, View, Canvas);
}

void FSpringMagicEditMode::DrawTextItem(const FText& Text, FCanvas* Canvas, float X, float& Y, float FontHeight)
{
	FCanvasTextItem TextItem(FVector2D::ZeroVector, Text, GEngine->GetSmallFont(), FLinearColor::White);
	TextItem.EnableShadow(FLinearColor::Black);
	Canvas->DrawItem(TextItem, X, Y);
	Y -= (3 + FontHeight);
}

void FSpringMagicEditMode::Draw3DTextItem(const FText& Text, FCanvas* Canvas, const FSceneView* View,
                                            const FViewport* Viewport, FVector Location)
{
	const int32 HalfX = Viewport->GetSizeXY().X / 2 / Canvas->GetDPIScale();
	const int32 HalfY = Viewport->GetSizeXY().Y / 2 / Canvas->GetDPIScale();

	const FPlane proj = View->Project(Location);
	if (proj.W > 0.f)
	{
		const int32 XPos = HalfX + (HalfX * proj.X);
		const int32 YPos = HalfY + (HalfY * (proj.Y * -1));
		FCanvasTextItem TextItem(FVector2D(XPos, YPos), Text, GEngine->GetSmallFont(), FLinearColor::White);
		TextItem.EnableShadow(FLinearColor::Black);
		Canvas->DrawItem(TextItem);
	}
}

#undef LOCTEXT_NAMESPACE
