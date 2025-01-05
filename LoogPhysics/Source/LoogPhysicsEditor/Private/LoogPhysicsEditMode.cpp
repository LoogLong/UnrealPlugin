// Copyright LoogLong. All Rights Reserved.
#include "LoogPhysicsEditMode.h"
#include "CanvasItem.h"
#include "CanvasTypes.h"
#include "EditorViewportClient.h"
#include "IPersonaPreviewScene.h"
#include "SceneManagement.h"
#include "Animation/DebugSkelMeshComponent.h"
#include "Materials/MaterialInstanceDynamic.h"

#define LOCTEXT_NAMESPACE "LoogPhysicsEditMode"


FEditorModeID FLoogPhysicsEditMode::ModeName = "AnimGraph.SkeletalControl.LoogPhysics";
FLoogPhysicsEditMode::FLoogPhysicsEditMode()
	: RuntimeNode(nullptr)
	  , GraphNode(nullptr)
{
}

void FLoogPhysicsEditMode::EnterMode(UAnimGraphNode_Base* InEditorNode, FAnimNode_Base* InRuntimeNode)
{
	RuntimeNode = static_cast<FAnimNode_LoogPhysics*>(InRuntimeNode);
	GraphNode = CastChecked<UAnimGraphNode_LoogPhysics>(InEditorNode);

	FAnimNodeEditMode::EnterMode(InEditorNode, InRuntimeNode);
}

void FLoogPhysicsEditMode::ExitMode()
{
	GraphNode = nullptr;
	RuntimeNode = nullptr;

	FAnimNodeEditMode::ExitMode();
}

void FLoogPhysicsEditMode::Render(const FSceneView* View, FViewport* Viewport, FPrimitiveDrawInterface* PDI)
{
	const USkeletalMeshComponent* SkelMeshComp = GetAnimPreviewScene().GetPreviewMeshComponent();

	if (RuntimeNode != nullptr)
	{
		const auto& RuntimeParticles = RuntimeNode->Particles;
		if (RuntimeParticles.Num() == 0)
		{
			return;
		}
		if (GraphNode->bDrawParticles)
		{
			int32 ParticleCount = RuntimeParticles.Num();
			for (int32 ParticleIndex = 0; ParticleIndex < ParticleCount; ++ParticleIndex)
			{
				const auto& RuntimeParticle = RuntimeParticles[ParticleIndex];

				DrawWireSphere(PDI, RuntimeParticle.Position, FLinearColor::Blue, RuntimeParticle.Thickness, 16, SDPG_Foreground);
			}
		}
		if (GraphNode->bDrawConstraints)
		{
			for (const auto& Section : RuntimeNode->ClothSections)
			{
				for (auto& Constraint : Section.LocalVerticalConstraints)
				{
					const auto& AParticle = RuntimeParticles[Constraint.ParticleAIndex];
					const auto& BParticle = RuntimeParticles[Constraint.ParticleBIndex];
					PDI->DrawLine(AParticle.Position, BParticle.Position, FLinearColor::White, SDPG_Foreground);
				}
				for (auto& Constraint : Section.GlobalStructureConstraints)
				{
					const auto& AParticle = RuntimeParticles[Constraint.ParticleAIndex];
					const auto& BParticle = RuntimeParticles[Constraint.ParticleBIndex];
					PDI->DrawLine(AParticle.Position, BParticle.Position, FLinearColor::White, SDPG_Foreground);
				}
				for (auto& Constraint : Section.LocalHorizontalConstraints)
				{
					const auto& AParticle = RuntimeParticles[Constraint.ParticleAIndex];
					const auto& BParticle = RuntimeParticles[Constraint.ParticleBIndex];
					PDI->DrawLine(AParticle.Position, BParticle.Position, FLinearColor::White, SDPG_Foreground);
				}
				for (auto& Constraint : Section.BendingStructureConstraints)
				{
					const auto& AParticle = RuntimeParticles[Constraint.ParticleAIndex];
					const auto& BParticle = RuntimeParticles[Constraint.ParticleBIndex];
					PDI->DrawLine(AParticle.Position, BParticle.Position, FLinearColor::White, SDPG_Foreground);
				}
			}
		}
		if (GraphNode->bDrawColliders)
		{
			const FMaterialRenderProxy* MaterialRenderProxy = GEngine->ConstraintLimitMaterialPrismatic->GetRenderProxy();
			int32                       ColliderCount       = RuntimeNode->SetupColliders.Num();
			for (int32 Idx = 0; Idx < ColliderCount; ++Idx)
			{
				auto& SetupCollider   = RuntimeNode->SetupColliders[Idx];
				auto& RuntimeCollider = RuntimeNode->RuntimeColliders[Idx];

				SetupCollider.DrawElemSolid(PDI, RuntimeCollider.TransformSimSpace, 1.f, MaterialRenderProxy);
				SetupCollider.DrawElemWire(PDI, RuntimeCollider.TransformSimSpace, 1.f, FColor::Black);
			}
		}
		if (GraphNode->bDrawClothColliders)
		{
			const FMaterialRenderProxy* MaterialRenderProxy = GEngine->ConstraintLimitMaterialPrismatic->GetRenderProxy();
			for (const auto& Section : RuntimeNode->ClothSections)
			{
				for (auto& Constraint : Section.Colliders)
				{
					const auto& AParticle = RuntimeParticles[Constraint.Particle0Index];
					const auto& BParticle = RuntimeParticles[Constraint.Particle1Index];
					PDI->DrawLine(AParticle.Position, BParticle.Position, FLinearColor::Black, SDPG_Foreground);
					DrawCylinder(PDI, AParticle.Position, BParticle.Position, AParticle.Thickness, 25, MaterialRenderProxy, SDPG_World);
				}
			}
		}
	}
	FAnimNodeEditMode::Render(View, Viewport, PDI);
}

void FLoogPhysicsEditMode::DrawHUD(FEditorViewportClient* ViewportClient, FViewport* Viewport, const FSceneView* View,
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

void FLoogPhysicsEditMode::DrawTextItem(const FText& Text, FCanvas* Canvas, float X, float& Y, float FontHeight)
{
	FCanvasTextItem TextItem(FVector2D::ZeroVector, Text, GEngine->GetSmallFont(), FLinearColor::White);
	TextItem.EnableShadow(FLinearColor::Black);
	Canvas->DrawItem(TextItem, X, Y);
	Y -= (3 + FontHeight);
}

void FLoogPhysicsEditMode::Draw3DTextItem(const FText& Text, FCanvas* Canvas, const FSceneView* View,
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
