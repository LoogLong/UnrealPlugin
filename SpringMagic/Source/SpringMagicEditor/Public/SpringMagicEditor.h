// Copyright LoogLong. All Rights Reserved.
#pragma once

#include "CoreMinimal.h"

class FSpringMagicEditorModule : public IModuleInterface
{
public:
	/** IModuleInterface implementation */
	virtual void StartupModule() override;
	virtual void ShutdownModule() override;
};
