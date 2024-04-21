﻿using Godot;

/// <summary>
///   Direct starter for aware (basically <see cref="MainGameState.MulticellularStage"/> with aware species)
/// </summary>
public partial class AwareStageStarter : ComplexStageStarterBase
{
    protected override MainGameState SimplyLoadableGameState => MainGameState.MulticellularStage;

    protected override void CustomizeLoadedScene(Node scene)
    {
        var stage = (MulticellularStage)scene;

        // Setup a new game where the player is specifically an aware creature already
        var game = GameProperties.StartNewAwareStageGame(new WorldGenerationSettings());

        stage.CurrentGame = game;
    }
}
