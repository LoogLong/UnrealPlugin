// Copyright LoogLong. All Rights Reserved.
#pragma once
#include "AnimNode_LoogPhysics.h"
#include "AnimGraphNode_SkeletalControlBase.h"
#include "EdGraph/EdGraphNodeUtils.h"

#include "AnimGraphNode_LoogPhysics.generated.h"


USTRUCT()
struct FLoogPhysicsBoneChain
{
	GENERATED_BODY()

	UPROPERTY(EditAnywhere)
	FBoneReference RootBone;

	UPROPERTY(EditAnywhere)
	FBoneReference EndBone;

	UPROPERTY(EditAnywhere)
	float VirtualBoneLength = 10.0f;

	UPROPERTY(EditAnywhere)
	float Mass = 1.0f;

	UPROPERTY(EditAnywhere)
	float Thickness = 1.0f;
};


USTRUCT()
struct FLoogPhysicsBoneSection
{
	GENERATED_BODY()

	UPROPERTY(EditAnywhere)
	TArray<FLoogPhysicsBoneChain> BoneChains;

	UPROPERTY(EditAnywhere)
	bool bNeedConnectChains = true;

	UPROPERTY(EditAnywhere)
	bool bChainLoop = true;

	UPROPERTY(EditAnywhere)
	bool bConnectChainsCollision = true;
};

UCLASS()
class UAnimGraphNode_LoogPhysics : public UAnimGraphNode_SkeletalControlBase
{
	GENERATED_UCLASS_BODY()
	UPROPERTY(EditAnywhere, Category = Settings)
	FAnimNode_LoogPhysics Node;

	virtual FText GetNodeTitle(ENodeTitleType::Type TitleType) const override;

	// UObject interface
	virtual void PostEditChangeProperty(struct FPropertyChangedEvent& PropertyChangedEvent) override;

protected:
	// UAnimGraphNode_Base interface
	virtual FEditorModeID GetEditorMode() const override;
	// virtual void ValidateAnimNodePostCompile(FCompilerResultsLog& MessageLog,
	//                                          UAnimBlueprintGeneratedClass* CompiledClass,
	//                                          int32 CompiledNodeIndex) override;
	virtual void CopyNodeDataToPreviewNode(FAnimNode_Base* AnimNode) override;
	virtual void CustomizeDetails(IDetailLayoutBuilder& DetailBuilder) override;
	// End of UAnimGraphNode_Base interface

	//virtual FText GetControllerDescription() const override;
	virtual FText GetControllerDescription() const override;
	virtual const FAnimNode_SkeletalControlBase* GetNode() const override { return &Node; }
	// End of UAnimGraphNode_SkeletalControlBase interface

	// UObject interface
	virtual void Serialize(FArchive& Ar) override;

	// End of UObject interface

public:
	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	TObjectPtr<USkeletalMesh> SkeletalMesh;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	TArray<FLoogPhysicsBoneSection> BoneSections;

	// UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	// TObjectPtr<UPhysicsAsset> PhysicsAsset;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	float ShrinkCompliance = 0.01f;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	float StretchCompliance = 0.01f;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	bool bBendingVertical = true;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	bool bBendingHorizontal = true;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	FString PropertyName;

	UPROPERTY(EditAnywhere, Category = "Loog Physics Tools")
	TArray<float> ValueToSet;
private:
	/** Constructing FText strings can be costly, so we cache the node's title */
	FNodeTitleTextTable CachedNodeTitles;

	void CreateParticlesAndConstraints();

	void CreateShearConstraints();

	void CreateBendingConstraints();
	void RemoveConstraints();
	void ModifyParticleProperty();

	bool IsChainValid();
};
