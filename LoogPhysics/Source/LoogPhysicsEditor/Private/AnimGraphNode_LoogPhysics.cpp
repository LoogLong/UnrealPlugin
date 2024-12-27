// Copyright LoogLong. All Rights Reserved.
#include "AnimGraphNode_LoogPhysics.h"

#include "DetailCategoryBuilder.h"
#include "DetailLayoutBuilder.h"
#include "DetailWidgetRow.h"
#include "LoogPhysicsEditMode.h"
#include "Kismet2/CompilerResultsLog.h"
#include "Widgets/Layout/SUniformGridPanel.h"

#define LOCTEXT_NAMESPACE "LoogPhysics"

// ----------------------------------------------------------------------------
UAnimGraphNode_LoogPhysics::UAnimGraphNode_LoogPhysics(const FObjectInitializer& ObjectInitializer)
	: Super(ObjectInitializer)
{
}

FText UAnimGraphNode_LoogPhysics::GetControllerDescription() const
{
	return LOCTEXT("LoogPhysics", "Loog Physics");
}


// ----------------------------------------------------------------------------
FText UAnimGraphNode_LoogPhysics::GetNodeTitle(ENodeTitleType::Type TitleType) const
{
	return LOCTEXT("NodeTitle", "Loog Physics");
}

void UAnimGraphNode_LoogPhysics::PostEditChangeProperty(struct FPropertyChangedEvent& PropertyChangedEvent)
{
	Super::PostEditChangeProperty(PropertyChangedEvent);

}

FEditorModeID UAnimGraphNode_LoogPhysics::GetEditorMode() const
{
	return FLoogPhysicsEditMode::ModeName;
}

void UAnimGraphNode_LoogPhysics::CopyNodeDataToPreviewNode(FAnimNode_Base* AnimNode)
{
	FAnimNode_LoogPhysics* LoogPhysics = reinterpret_cast<FAnimNode_LoogPhysics*>(AnimNode);
	
	Node.ClothSections[0].Constraints;
	for (int i = 0; i < Node.ClothSections.Num(); ++i)
	{
		LoogPhysics->ClothSections[i].Constraints = Node.ClothSections[i].Constraints;
	}
}

void UAnimGraphNode_LoogPhysics::CustomizeDetails(IDetailLayoutBuilder& DetailBuilder)
{
	IDetailCategoryBuilder& ViewportCategory = DetailBuilder.EditCategory(TEXT("Loog Physics Tools"));
	FDetailWidgetRow& WidgetRow = ViewportCategory.AddCustomRow(LOCTEXT("LoogPhysics", "LoogPhysicsTools"));

	WidgetRow
		[
			SNew(SUniformGridPanel)
				.SlotPadding(FMargin(2, 0, 2, 0))
				// + SUniformGridPanel::Slot(0, 0)
				// [
				// 	SNew(SButton)
				// 		.HAlign(HAlign_Center)
				// 		.VAlign(VAlign_Center)
				// 		.OnClicked_Lambda([this]()
				// 			{
				// 				this->CreateStructureConstraints();
				// 				return FReply::Handled();
				// 			})
				// 		.Content()
				// 		[
				// 			SNew(STextBlock)
				// 				.Text(FText::FromString(TEXT("CreateStructureConstraints")))
				// 		]
				// ]
				+ SUniformGridPanel::Slot(0, 0)
				[
					SNew(SButton)
						.HAlign(HAlign_Center)
						.VAlign(VAlign_Center)
						.OnClicked_Lambda([this]()
							{
								this->CreateBendingConstraints();
								return FReply::Handled();
							})
						.Content()
						[
							SNew(STextBlock)
								.Text(FText::FromString(TEXT("CreateBendingConstraints")))
						]
				]
				+ SUniformGridPanel::Slot(0, 1)
				[
					SNew(SButton)
						.HAlign(HAlign_Center)
						.VAlign(VAlign_Center)
						.OnClicked_Lambda([this]()
							{
								this->CreateShearConstraints();
								return FReply::Handled();
							})
						.Content()
						[
							SNew(STextBlock)
								.Text(FText::FromString(TEXT("CreateShearingConstraints")))
						]
				]
				+ SUniformGridPanel::Slot(0, 2)
				[
					SNew(SButton)
						.HAlign(HAlign_Center)
						.VAlign(VAlign_Center)
						.OnClicked_Lambda([this]()
							{
								this->RemoveConstraints();
								return FReply::Handled();
							})
						.Content()
						[
							SNew(STextBlock)
								.Text(FText::FromString(TEXT("RemoveConstraints")))
						]
				]
				+ SUniformGridPanel::Slot(0, 3)
				[
					SNew(SButton)
						.HAlign(HAlign_Center)
						.VAlign(VAlign_Center)
						.OnClicked_Lambda([this]()
							{
								this->ModifyParticleProperty();
								return FReply::Handled();
							})
						.Content()
						[
							SNew(STextBlock)
								.Text(FText::FromString(TEXT("ModifyParticleProperty")))
						]
				]
		];
	Super::CustomizeDetails(DetailBuilder);
}


void UAnimGraphNode_LoogPhysics::Serialize(FArchive& Ar)
{
	Super::Serialize(Ar);
}
void UAnimGraphNode_LoogPhysics::CreateStructureConstraints()
{
	for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
	{



	}
}

void UAnimGraphNode_LoogPhysics::CreateShearConstraints()
{
	if (!IsChainValid())
	{
		FString ErrorString = FString::Printf(TEXT("Chain Length Must The Same!"));
		FMessageDialog::Open(EAppMsgType::Ok, FText::FromString(ErrorString));
		return;
	}

	for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
	{
		int32 ChainNum = ClothSection.Chains.Num();

		for (int32 ChainIndex = 0; ChainIndex < ChainNum; ++ChainIndex)
		{
			auto& Chain = ClothSection.Chains[ChainIndex];
			int32 ChainLength = Chain.ParticleIndices.Num();

			int32 NextChainIndex = ChainIndex + 1;
			if (NextChainIndex < ChainNum || ClothSection.bChainLoop)
			{
				if (NextChainIndex >= ChainNum)
				{
					NextChainIndex = 0;
				}
				auto& NextChain = ClothSection.Chains[NextChainIndex];
				for (int32 Idx = 1; Idx < ChainLength; ++Idx)
				{
					{// left top - right bottom
						int32 AIndex = Chain.ParticleIndices[Idx - 1];
						int32 BIndex = NextChain.ParticleIndices[Idx];

						auto& Constraint = ClothSection.Constraints.AddDefaulted_GetRef();
						Constraint.ParticleAIndex = AIndex;
						Constraint.ParticleBIndex = BIndex;
						Constraint.ShrinkCompliance = ShrinkCompliance;
						Constraint.StretchCompliance = StretchCompliance;
					}
					{ // left bottom - right top
						int32 AIndex = Chain.ParticleIndices[Idx];
						int32 BIndex = NextChain.ParticleIndices[Idx - 1];

						auto& Constraint = ClothSection.Constraints.AddDefaulted_GetRef();
						Constraint.ParticleAIndex = AIndex;
						Constraint.ParticleBIndex = BIndex;
						Constraint.ShrinkCompliance = ShrinkCompliance;
						Constraint.StretchCompliance = StretchCompliance;
					}
				}
			}
		}
	}
}

void UAnimGraphNode_LoogPhysics::CreateBendingConstraints()
{
	if (!IsChainValid())
	{
		FString ErrorString = FString::Printf(TEXT("Chain Length Must The Same!"));
		FMessageDialog::Open(EAppMsgType::Ok, FText::FromString(ErrorString));
		return;
	}
	for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
	{
		int32 ChainNum = ClothSection.Chains.Num();

		for (int32 ChainIndex = 0; ChainIndex < ChainNum; ++ChainIndex)
		{
			auto& Chain = ClothSection.Chains[ChainIndex];
			int32 ChainLength = Chain.ParticleIndices.Num();

			if (bBendingVertical)
			{
				for (int32 Idx = 2; Idx < ChainLength; ++Idx)
				{
					int32 AIndex = Chain.ParticleIndices[Idx - 2];
					int32 BIndex = Chain.ParticleIndices[Idx];

					auto& Constraint = ClothSection.Constraints.AddDefaulted_GetRef();
					Constraint.ParticleAIndex = AIndex;
					Constraint.ParticleBIndex = BIndex;
					Constraint.ShrinkCompliance = ShrinkCompliance;
					Constraint.StretchCompliance = StretchCompliance;
				}
			}
			
			if (bBendingHorizontal)
			{
				int32 NextChainIndex = ChainIndex + 2;
				if (NextChainIndex < ChainNum || ClothSection.bChainLoop)
				{
					if (NextChainIndex >= ChainNum)
					{
						NextChainIndex = FMath::Modulo(NextChainIndex, ChainNum);
					}
					auto& NextChain = ClothSection.Chains[NextChainIndex];
					for (int32 Idx = 1; Idx < ChainLength; ++Idx)
					{
						int32 AIndex = Chain.ParticleIndices[Idx];
						int32 BIndex = NextChain.ParticleIndices[Idx];

						auto& Constraint = ClothSection.Constraints.AddDefaulted_GetRef();
						Constraint.ParticleAIndex = AIndex;
						Constraint.ParticleBIndex = BIndex;
						Constraint.ShrinkCompliance = ShrinkCompliance;
						Constraint.StretchCompliance = StretchCompliance;
					}
				}
			}
			
		}
	}
}

void UAnimGraphNode_LoogPhysics::RemoveConstraints()
{
	for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
	{
		ClothSection.Constraints.Empty();
	}
}

void UAnimGraphNode_LoogPhysics::ModifyParticleProperty()
{
	if (PropertyName.IsEmpty())
	{
		return;
	}
	if (ValueToSet.IsEmpty())
	{
		return;
	}
	auto ParticleStruct = FLoogPhysicsParticle::StaticStruct();
	auto ParticleStructProperty = ParticleStruct->FindPropertyByName(FName(PropertyName));
	if (ParticleStructProperty == nullptr)
	{
		return;
	}
	bool bSingleValue = ValueToSet.Num() == 1;
	if (bSingleValue)
	{
		for (FLoogPhysicsParticle& Particle : Node.Particles)
		{
			ParticleStructProperty->SetValue_InContainer(&Particle, &ValueToSet[0]);
		}
	}
	else
	{
		// chain value
		if (IsChainValid())
		{
			bool bValueNumValid = Node.ClothSections[0].Chains[0].ParticleIndices.Num() == ValueToSet.Num();
			if (bValueNumValid)
			{
				for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
				{
					for (FLoogPhysicsChain& Chain : ClothSection.Chains)
					{
						for (int ValueIndex = 0; ValueIndex < Chain.ParticleIndices.Num(); ++ValueIndex)
						{
							int32& ParticleIndex = Chain.ParticleIndices[ValueIndex];
							FLoogPhysicsParticle& Particle = Node.Particles[ParticleIndex];
							ParticleStructProperty->SetValue_InContainer(&Particle, &ValueToSet[ValueIndex]);
						}
					}
				}
			}
		}
	}
}

bool UAnimGraphNode_LoogPhysics::IsChainValid()
{
	for (FLoogPhysicsClothSection& ClothSection : Node.ClothSections)
	{
		int32 ChainLength = ClothSection.Chains[0].ParticleIndices.Num();
		for (int32 i = 1; i < ClothSection.Chains.Num(); ++i)
		{
			if (ClothSection.Chains[i].ParticleIndices.Num() != ChainLength)
			{
				return false;
			}
		}
	}
	return true;
}

#undef LOCTEXT_NAMESPACE
