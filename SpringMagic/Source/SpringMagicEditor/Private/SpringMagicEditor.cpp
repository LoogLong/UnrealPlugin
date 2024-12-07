// Copyright LoogLong. All Rights Reserved.
#include "SpringMagicEditor.h"
#include "Modules/ModuleManager.h"
#include "Textures/SlateIcon.h"
#include "SpringMagicEditMode.h"

#define LOCTEXT_NAMESPACE "SpringMagicEditorModule"


void FSpringMagicEditorModule::StartupModule()
{
	FEditorModeRegistry::Get().RegisterMode<FSpringMagicEditMode>(FSpringMagicEditMode::ModeName, LOCTEXT("FSpringMagicEditMode", "Spring Magic"), FSlateIcon(), false);
}


void FSpringMagicEditorModule::ShutdownModule()
{
	FEditorModeRegistry::Get().UnregisterMode(FSpringMagicEditMode::ModeName);
}

#undef LOCTEXT_NAMESPACE

IMPLEMENT_MODULE(FSpringMagicEditorModule, SpringMagicEditor)
