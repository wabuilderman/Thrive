using System;
using Godot;
using Newtonsoft.Json;

/// <summary>
///   Pulse of particles representing the player cell's chemoreception ability
/// </summary>
[JSONAlwaysDynamicType]
[SceneLoadedClass("res://src/microbe_stage/ChemoreceptorPulse.tscn", UsesEarlyResolve = false)]
public class ChemoreceptorPulse : RigidBody, ITimedLife
{
    private Particles particles = null!;

    public float TimeToLiveRemaining { get; set; }

    [JsonProperty]
    private float? FadeTimeRemaining { get; set; }

    public EntityReference<IEntity> Emitter { get; set; } = new();

    public Color Colour { get; set; } = new();

    public void OnTimeOver()
    {
        if (FadeTimeRemaining == null)
            BeginDestroy();
    }

    public override void _Ready()
    {
        particles = GetNode<Particles>("Particles");
        particles.ProcessMaterial.Set("color", Colour);

        var emitterNode = Emitter.Value?.EntityNode;

        if (emitterNode != null)
            AddCollisionExceptionWith(emitterNode);
    }

    public override void _Process(float delta)
    {
        if (FadeTimeRemaining == null)
            return;

        FadeTimeRemaining -= delta;
        if (FadeTimeRemaining <= 0)
            Destroy();
    }

    /// <summary>
    ///   Stops particle emission and destroys the object after 5 seconds.
    /// </summary>
    private void BeginDestroy()
    {
        particles.Emitting = false;
        CollisionLayer = 0;
        CollisionMask = 0;
        LinearVelocity = Vector3.Zero;

        // Timer that delays despawn of projectiles
        FadeTimeRemaining = Constants.PROJECTILE_DESPAWN_DELAY;
    }

    private void Destroy()
    {
        this.DetachAndQueueFree();
    }
}