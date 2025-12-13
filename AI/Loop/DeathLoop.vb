Imports System.Collections.Generic
Imports System.Drawing
Imports System.Drawing.Drawing2D
Imports System.Numerics
Imports System.Windows.Forms

Public Class DeathLoopGame
    Inherits Form
    
    Private WithEvents gameTimer As New Timer()
    Private WithEvents physicsTimer As New Timer()
    Private renderTarget As Bitmap
    Private renderGraphics As Graphics
    
    ' Game State
    Private gameTime As Double = 0
    Private loopTime As Double = 0
    Private loopCount As Integer = 0
    Private gameState As GameState = GameState.Running
    Private playerHasDeathLoopPower As Boolean = True
    
    ' Physics World
    Private physicsWorld As New PhysicsWorld()
    Private ragdolls As New List(Of RagdollPhysics)()
    Private destructibleObjects As New List(Of DestructibleObject)()
    Private timeAnomalies As New List(Of TimeAnomaly)()
    Private gravityWells As New List(Of GravityWell)()
    
    ' Player
    Private player As New DeathLoopPlayer()
    Private playerInventory As New Inventory()
    Private playerMemory As New LoopMemory()
    
    ' Visionaries (Targets)
    Private visionaries As New List(Of Visionary)()
    
    ' Game Events
    Private worldEvents As New EventSystem()
    Private loopEvents As New LoopEventManager()
    
    ' Constants
    Private Const LOOP_DURATION As Double = 7200 ' 2 hours in seconds
    Private Const PHYSICS_TIME_STEP As Double = 1.0 / 120.0 ' 120Hz physics
    Private Const RENDER_TIME_STEP As Double = 1.0 / 60.0 ' 60Hz rendering
    
    Public Enum GameState
        Running
        Paused
        LoopReset
        PlayerDead
        VisionaryKilled
    End Enum
    
    ' ============================================
    ' PHYSICS SYSTEM
    ' ============================================
    
    Public Class PhysicsWorld
        Public Property Gravity As Vector2 = New Vector2(0, 98.1F) ' 10x Earth gravity for game feel
        Public Property AirDensity As Single = 1.225F ' kg/m³
        Public Property FrictionCoefficient As Single = 0.3F
        Public Property Restitution As Single = 0.2F ' Bounciness
        Public Property TimeDilation As Single = 1.0F
        
        Private collisionPairs As New List(Of CollisionPair)()
        Private spatialGrid As New SpatialHashGrid(50) ' 50px cell size
        
        Public Sub Update(deltaTime As Double)
            ' Apply temporal effects
            Dim effectiveDeltaTime As Single = deltaTime * TimeDilation
            
            ' Update spatial grid
            spatialGrid.Clear()
            For Each obj In physicsObjects
                spatialGrid.Insert(obj)
            Next
            
            ' Detect collisions
            DetectCollisions()
            
            ' Resolve collisions
            For Each pair In collisionPairs
                ResolveCollision(pair)
            Next
            
            ' Clear collision pairs for next frame
            collisionPairs.Clear()
        End Sub
        
        Private Sub DetectCollisions()
            ' Broad phase using spatial hashing
            For Each cell In spatialGrid.GetNonEmptyCells()
                Dim objectsInCell = spatialGrid.GetObjectsInCell(cell)
                
                For i As Integer = 0 To objectsInCell.Count - 2
                    For j As Integer = i + 1 To objectsInCell.Count - 1
                        Dim objA = objectsInCell(i)
                        Dim objB = objectsInCell(j)
                        
                        If ShouldCollide(objA, objB) AndAlso
                           CheckCollision(objA, objB) Then
                            collisionPairs.Add(New CollisionPair(objA, objB))
                        End If
                    Next
                Next
            Next
        End Sub
        
        Private Function ShouldCollide(objA As PhysicsObject, objB As PhysicsObject) As Boolean
            ' Ignore collisions between objects with the same faction unless enabled
            If objA.Faction = objB.Faction AndAlso Not objA.CollideWithSameFaction Then
                Return False
            End If
            
            ' Check collision layers
            Return (objA.CollisionLayer And objB.CollisionMask) <> 0 AndAlso
                   (objB.CollisionLayer And objA.CollisionMask) <> 0
        End Function
        
        Private Function CheckCollision(objA As PhysicsObject, objB As PhysicsObject) As Boolean
            ' GJK (Gilbert-Johnson-Keerthi) algorithm for convex shapes
            Return GJKIntersection(objA, objB)
        End Function
        
        Private Sub ResolveCollision(pair As CollisionPair)
            Dim objA = pair.ObjectA
            Dim objB = pair.ObjectB
            
            ' Calculate collision normal and penetration depth
            Dim collisionInfo As CollisionInfo = CalculateCollisionInfo(objA, objB)
            
            If collisionInfo Is Nothing Then Return
            
            ' Separate objects
            Dim totalInverseMass As Single = objA.InverseMass + objB.InverseMass
            If totalInverseMass > 0 Then
                Dim separation = collisionInfo.Penetration / totalInverseMass
                objA.Position -= collisionInfo.Normal * separation * objA.InverseMass
                objB.Position += collisionInfo.Normal * separation * objB.InverseMass
            End If
            
            ' Calculate relative velocity
            Dim relativeVelocity = objB.Velocity - objA.Velocity
            Dim velocityAlongNormal = Vector2.Dot(relativeVelocity, collisionInfo.Normal)
            
            ' Do not resolve if velocities are separating
            If velocityAlongNormal > 0 Then Return
            
            ' Calculate restitution (bounciness)
            Dim e = Math.Min(objA.Restitution, objB.Restitution)
            
            ' Calculate impulse scalar
            Dim j = -(1 + e) * velocityAlongNormal
            j /= totalInverseMass
            
            ' Apply impulse
            Dim impulse = collisionInfo.Normal * j
            objA.Velocity -= impulse * objA.InverseMass
            objB.Velocity += impulse * objB.InverseMass
            
            ' Apply friction
            ApplyFriction(objA, objB, collisionInfo, impulse)
            
            ' Fire collision events
            objA.OnCollision(objB, collisionInfo)
            objB.OnCollision(objA, collisionInfo)
        End Sub
        
        Private Sub ApplyFriction(objA As PhysicsObject, objB As PhysicsObject, 
                                 collisionInfo As CollisionInfo, impulse As Vector2)
            ' Calculate tangent vector
            Dim relativeVelocity = objB.Velocity - objA.Velocity
            Dim tangent = relativeVelocity - collisionInfo.Normal * 
                         Vector2.Dot(relativeVelocity, collisionInfo.Normal)
            
            If tangent.LengthSquared() > 0.0001F Then
                tangent = Vector2.Normalize(tangent)
                
                ' Magnitude to apply along friction vector
                Dim jt = -Vector2.Dot(relativeVelocity, tangent)
                jt /= (objA.InverseMass + objB.InverseMass)
                
                ' Coulomb's law: Friction <= μ * N
                Dim mu = Math.Sqrt(objA.Friction * objA.Friction + 
                                  objB.Friction * objB.Friction)
                Dim frictionImpulse As Vector2
                
                If Math.Abs(jt) < impulse.Length() * mu Then
                    frictionImpulse = tangent * jt
                Else
                    frictionImpulse = tangent * -impulse.Length() * mu
                End If
                
                ' Apply friction impulse
                objA.Velocity -= frictionImpulse * objA.InverseMass
                objB.Velocity += frictionImpulse * objB.InverseMass
            End If
        End Sub
    End Class
    
    ' ============================================
    ' PLAYER WITH TIME LOOP ABILITIES
    ' ============================================
    
    Public Class DeathLoopPlayer
        Inherits PhysicsObject
        
        Public Property Health As Single = 100.0F
        Public Property MaxHealth As Single = 100.0F
        Public Property Energy As Single = 100.0F
        Public Property MaxEnergy As Single = 100.0F
        Public Property IsInvulnerable As Boolean = False
        Public Property InvulnerabilityTimer As Double = 0
        
        ' Time Loop Powers
        Public Property CanRewindTime As Boolean = False
        Public Property CanSlowTime As Boolean = False
        Public Property CanFreezeTime As Boolean = False
        Public Property CanCreateEcho As Boolean = False
        Public Property CanPhaseShift As Boolean = False
        
        ' Movement
        Public Property MoveSpeed As Single = 5.0F
        Public Property SprintMultiplier As Single = 1.5F
        Public Property JumpForce As Single = 500.0F
        Public Property IsGrounded As Boolean = False
        Public Property IsSprinting As Boolean = False
        Public Property IsCrouching As Boolean = False
        
        ' Combat
        Public Property CurrentWeapon As Weapon = Nothing
        Public Property MeleeDamage As Single = 25.0F
        Public Property CriticalChance As Single = 0.1F
        Public Property CriticalMultiplier As Single = 2.0F
        
        ' Loop Memory
        Public Property RememberedPaths As New List(Of Vector2)()
        Public Property LoopDeaths As Integer = 0
        Public Property LoopKills As Integer = 0
        Public Property DiscoveredSecrets As New List(Of String)()
        
        ' Temporal Echo (record/replay)
        Private echoRecording As New List(Of PlayerFrame)()
        Private echoPlayingBack As Boolean = False
        Private echoPlaybackIndex As Integer = 0
        
        Public Sub Update(deltaTime As Double, input As PlayerInput)
            ' Handle time powers
            UpdateTimePowers(deltaTime, input)
            
            ' Handle movement
            UpdateMovement(deltaTime, input)
            
            ' Handle combat
            UpdateCombat(deltaTime, input)
            
            ' Update physics
            MyBase.Update(deltaTime)
            
            ' Record for echo if active
            If CanCreateEcho AndAlso input.IsRecordingEcho Then
                RecordEchoFrame()
            End If
        End Sub
        
        Private Sub UpdateTimePowers(deltaTime As Double, input As PlayerInput)
            ' Time Rewind
            If CanRewindTime AndAlso input.IsRewindingTime Then
                If Energy > 0 Then
                    Energy -= 10.0F * deltaTime
                    ' Apply local time rewind effect
                    PhysicsWorld.TimeDilation = -0.5F ' Reverse time at half speed
                End If
            ElseIf CanSlowTime AndAlso input.IsSlowingTime Then
                If Energy > 0 Then
                    Energy -= 5.0F * deltaTime
                    PhysicsWorld.TimeDilation = 0.25F ' Quarter speed
                End If
            ElseIf CanFreezeTime AndAlso input.IsFreezingTime Then
                If Energy > 0 Then
                    Energy -= 20.0F * deltaTime
                    PhysicsWorld.TimeDilation = 0.001F ' Nearly stopped
                End If
            Else
                ' Normal time flow
                PhysicsWorld.TimeDilation = 1.0F
            End If
            
            ' Phase Shift
            If CanPhaseShift AndAlso input.IsPhasing Then
                If Energy > 0 Then
                    Energy -= 15.0F * deltaTime
                    Me.CollisionLayer = CollisionLayer.PhaseShifted
                    Me.Opacity = 0.5F
                End If
            Else
                Me.CollisionLayer = CollisionLayer.Player
                Me.Opacity = 1.0F
            End If
            
            ' Energy regeneration
            If Energy < MaxEnergy Then
                Energy += 2.0F * deltaTime
            End If
        End Sub
        
        Private Sub UpdateMovement(deltaTime As Double, input As PlayerInput)
            Dim moveDirection As Vector2 = Vector2.Zero
            
            ' Calculate movement direction
            If input.MoveUp Then moveDirection.Y -= 1
            If input.MoveDown Then moveDirection.Y += 1
            If input.MoveLeft Then moveDirection.X -= 1
            If input.MoveRight Then moveDirection.X += 1
            
            ' Normalize diagonal movement
            If moveDirection.LengthSquared() > 1 Then
                moveDirection = Vector2.Normalize(moveDirection)
            End If
            
            ' Apply sprint
            Dim speedMultiplier As Single = MoveSpeed
            If input.IsSprinting AndAlso Energy > 0 Then
                speedMultiplier *= SprintMultiplier
                IsSprinting = True
                Energy -= 5.0F * deltaTime
            Else
                IsSprinting = False
            End If
            
            ' Apply crouch
            If input.IsCrouching Then
                speedMultiplier *= 0.5F
                IsCrouching = True
            Else
                IsCrouching = False
            End If
            
            ' Apply movement
            If moveDirection.LengthSquared() > 0 Then
                Dim velocityChange = moveDirection * speedMultiplier
                Me.Velocity = New Vector2(velocityChange.X, Me.Velocity.Y)
                
                ' Remember path for loop memory
                RememberedPaths.Add(Me.Position)
                If RememberedPaths.Count > 1000 Then
                    RememberedPaths.RemoveAt(0)
                End If
            End If
            
            ' Jump
            If input.IsJumping AndAlso IsGrounded Then
                Me.Velocity.Y -= JumpForce
                IsGrounded = False
            End If
            
            ' Apply gravity
            If Not IsGrounded Then
                Me.Velocity.Y += PhysicsWorld.Gravity.Y * deltaTime
            End If
        End Sub
        
        Private Sub UpdateCombat(deltaTime As Double, input As PlayerInput)
            ' Melee attack
            If input.IsMeleeAttacking Then
                PerformMeleeAttack()
            End If
            
            ' Ranged attack
            If input.IsShooting AndAlso CurrentWeapon IsNot Nothing Then
                If CurrentWeapon.CanFire() Then
                    CurrentWeapon.Fire(Me.Position, input.AimDirection)
                End If
            End If
            
            ' Weapon switching
            If input.SwitchWeaponForward Then
                ' Cycle to next weapon
            ElseIf input.SwitchWeaponBackward Then
                ' Cycle to previous weapon
            End If
        End Sub
        
        Private Sub PerformMeleeAttack()
            ' Create melee hitbox
            Dim hitbox As New PhysicsObject()
            hitbox.Position = Me.Position + New Vector2(50, 0) ' Front of player
            hitbox.Size = New Vector2(100, 50)
            hitbox.CollisionLayer = CollisionLayer.PlayerAttack
            
            ' Check for collisions
            ' ... collision detection logic
            
            ' Apply damage to hit targets
            ' ... damage application logic
        End Sub
        
        Private Sub RecordEchoFrame()
            Dim frame As New PlayerFrame With {
                .Position = Me.Position,
                .Velocity = Me.Velocity,
                .Rotation = Me.Rotation,
                .Action = CurrentAction,
                .Timestamp = gameTime
            }
            echoRecording.Add(frame)
            
            ' Limit recording length
            If echoRecording.Count > 600 Then ' 10 seconds at 60fps
                echoRecording.RemoveAt(0)
            End If
        End Sub
        
        Public Sub PlayEcho()
            If echoRecording.Count > 0 AndAlso Not echoPlayingBack Then
                echoPlayingBack = True
                echoPlaybackIndex = 0
            End If
        End Sub
        
        Public Sub TakeDamage(damage As Single, damageType As DamageType, attacker As GameObject)
            If IsInvulnerable Then Return
            
            ' Calculate final damage with resistances
            Dim finalDamage = CalculateDamage(damage, damageType)
            Health -= finalDamage
            
            ' Apply knockback
            If attacker IsNot Nothing Then
                Dim direction = Vector2.Normalize(Me.Position - attacker.Position)
                Me.Velocity += direction * finalDamage * 10.0F
            End If
            
            ' Invulnerability frames
            If finalDamage > 0 Then
                IsInvulnerable = True
                InvulnerabilityTimer = 1.0 ' 1 second invulnerability
            End If
            
            ' Check for death
            If Health <= 0 Then
                OnDeath()
            End If
        End Sub
        
        Private Function CalculateDamage(baseDamage As Single, damageType As DamageType) As Single
            ' Apply resistances based on damage type
            Dim multiplier As Single = 1.0F
            
            ' TODO: Implement damage type resistances
            
            Return baseDamage * multiplier
        End Function
        
        Private Sub OnDeath()
            ' Record death in loop memory
            LoopDeaths += 1
            
            ' Create death ragdoll
            CreateRagdoll()
            
            ' Trigger loop reset if power is active
            If playerHasDeathLoopPower Then
                ' Reset to start of loop with memories intact
                ResetLoop()
            Else
                ' Game over
                gameState = GameState.PlayerDead
            End If
        End Sub
    End Class
    
    ' ============================================
    ' VISIONARY (TARGET) AI
    ' ============================================
    
    Public Class Visionary
        Inherits PhysicsObject
        
        Public Property Name As String
        Public Property Health As Single = 500.0F
        Public Property MaxHealth As Single = 500.0F
        Public Property Awareness As Single = 0.0F ' 0-100
        Public Property Aggression As Single = 50.0F ' 0-100
        Public Property Routine As New DailyRoutine()
        Public Property IsAlive As Boolean = True
        Public Property KillMethod As KillMethod = KillMethod.None
        Public Property UniqueAbility As VisionaryAbility
        
        ' AI States
        Public Property CurrentState As AIState = AIState.Patrolling
        Public Property TargetPosition As Vector2 = Vector2.Zero
        Public Property LastKnownPlayerPosition As Vector2 = Vector2.Zero
        Public Property SuspicionTimer As Double = 0
        
        ' Weapons and Abilities
        Public Property Weapons As New List(Of Weapon)()
        Public Property CurrentWeapon As Weapon = Nothing
        Public Property SpecialAbilityCooldown As Double = 0
        
        ' Loop Memory
        Public Property RemembersPlayer As Boolean = False
        Public Property LoopEncounters As Integer = 0
        
        Public Sub Update(deltaTime As Double, playerPosition As Vector2)
            ' Update AI state machine
            UpdateAIState(deltaTime, playerPosition)
            
            ' Update routine based on loop time
            Routine.Update(loopTime)
            
            ' Update awareness
            UpdateAwareness(deltaTime, playerPosition)
            
            ' Update physics
            MyBase.Update(deltaTime)
        End Sub
        
        Private Sub UpdateAIState(deltaTime As Double, playerPosition As Vector2)
            Select Case CurrentState
                Case AIState.Patrolling
                    UpdatePatrolling(deltaTime)
                    
                Case AIState.Investigating
                    UpdateInvestigating(deltaTime)
                    
                Case AIState.Combat
                    UpdateCombat(deltaTime, playerPosition)
                    
                Case AIState.Fleeing
                    UpdateFleeing(deltaTime, playerPosition)
                    
                Case AIState.UsingAbility
                    UpdateAbilityUsage(deltaTime)
            End Select
        End Sub
        
        Private Sub UpdatePatrolling(deltaTime As Double)
            ' Follow daily routine
            Dim routinePosition = Routine.GetCurrentPosition(loopTime)
            If routinePosition.HasValue Then
                TargetPosition = routinePosition.Value
                MoveToTarget(deltaTime)
            End If
            
            ' Check for player detection
            If Awareness > 50 Then
                CurrentState = AIState.Investigating
                SuspicionTimer = 5.0 ' 5 seconds to investigate
            End If
        End Sub
        
        Private Sub UpdateInvestigating(deltaTime As Double)
            ' Move to last known player position
            MoveToTarget(deltaTime)
            
            ' Look for player
            SuspicionTimer -= deltaTime
            
            If SuspicionTimer <= 0 Then
                If Awareness > 75 Then
                    CurrentState = AIState.Combat
                Else
                    CurrentState = AIState.Patrolling
                    Awareness = Math.Max(0, Awareness - 25)
                End If
            End If
        End Sub
        
        Private Sub UpdateCombat(deltaTime As Double, playerPosition As Vector2)
            ' Engage player
            LastKnownPlayerPosition = playerPosition
            TargetPosition = playerPosition
            
            ' Move toward player or take cover
            If Vector2.Distance(Position, playerPosition) > 200 Then
                MoveToTarget(deltaTime)
            Else
                ' Take cover and shoot
                If CurrentWeapon IsNot Nothing AndAlso CurrentWeapon.CanFire() Then
                    CurrentWeapon.Fire(Position, Vector2.Normalize(playerPosition - Position))
                End If
            End If
            
            ' Check health for fleeing
            If Health < MaxHealth * 0.3 AndAlso Aggression < 70 Then
                CurrentState = AIState.Fleeing
            End If
            
            ' Use special ability
            If UniqueAbility IsNot Nothing AndAlso SpecialAbilityCooldown <= 0 Then
                CurrentState = AIState.UsingAbility
            End If
        End Sub
        
        Private Sub UpdateFleeing(deltaTime As Double, playerPosition As Vector2)
            ' Move away from player
            Dim fleeDirection = Vector2.Normalize(Position - playerPosition)
            TargetPosition = Position + fleeDirection * 1000
            MoveToTarget(deltaTime)
            
            ' Regain courage
            If Vector2.Distance(Position, playerPosition) > 500 Then
                Aggression += 10 * deltaTime
                If Aggression > 50 Then
                    CurrentState = AIState.Combat
                End If
            End If
        End Sub
        
        Private Sub UpdateAbilityUsage(deltaTime As Double)
            ' Use unique ability
            If UniqueAbility IsNot Nothing Then
                UniqueAbility.Activate(Me, player)
                SpecialAbilityCooldown = UniqueAbility.Cooldown
            End If
            
            ' Return to combat
            CurrentState = AIState.Combat
        End Sub
        
        Private Sub MoveToTarget(deltaTime As Double)
            Dim direction = Vector2.Normalize(TargetPosition - Position)
            If Single.IsNaN(direction.X) Then direction = Vector2.Zero
            
            Dim moveSpeed As Single = 2.0F
            If CurrentState = AIState.Fleeing Then
                moveSpeed = 3.0F
            End If
            
            Velocity = direction * moveSpeed
        End Sub
        
        Private Sub UpdateAwareness(deltaTime As Double, playerPosition As Vector2)
            Dim distanceToPlayer = Vector2.Distance(Position, playerPosition)
            
            ' Base awareness decay
            Awareness -= 5.0F * deltaTime
            
            ' Increase awareness if player is close
            If distanceToPlayer < 300 Then
                Dim visibility = CalculateVisibility(playerPosition)
                Awareness += visibility * 20.0F * deltaTime
                
                ' If player is very close or attacking, max awareness
                If distanceToPlayer < 100 OrElse player.CurrentWeapon.IsFiring Then
                    Awareness = Math.Min(100, Awareness + 50 * deltaTime)
                End If
            End If
            
            ' Cap awareness
            Awareness = Math.Max(0, Math.Min(100, Awareness))
            
            ' Remember player across loops after enough encounters
            If Awareness > 90 Then
                LoopEncounters += 1
                If LoopEncounters >= 3 Then
                    RemembersPlayer = True
                End If
            End If
        End Sub
        
        Private Function CalculateVisibility(playerPosition As Vector2) As Single
            ' Simple line-of-sight check
            ' In full implementation, use raycasting against obstacles
            Return 1.0F ' Placeholder
        End Function
        
        Public Sub TakeDamage(damage As Single, damageType As DamageType, attacker As DeathLoopPlayer)
            Health -= damage
            
            ' Increase awareness
            Awareness = Math.Min(100, Awareness + 30)
            
            ' Set combat state
            If CurrentState <> AIState.Combat Then
                CurrentState = AIState.Combat
                LastKnownPlayerPosition = attacker.Position
            End If
            
            ' Check for death
            If Health <= 0 Then
                OnDeath(attacker)
            End If
        End Sub
        
        Private Sub OnDeath(attacker As DeathLoopPlayer)
            IsAlive = False
            KillMethod = DetermineKillMethod(attacker)
            
            ' Record kill in player's memory
            attacker.LoopKills += 1
            
            ' Drop loot
            DropLoot()
            
            ' Create ragdoll with physics
            CreateRagdoll()
            
            ' Trigger world state change
            worldEvents.Trigger("VisionaryKilled", Name)
        End Sub
        
        Private Function DetermineKillMethod(attacker As DeathLoopPlayer) As KillMethod
            ' Determine special kill method based on circumstances
            If attacker.CurrentWeapon IsNot Nothing Then
                Return KillMethod.WeaponSpecific
            End If
            
            ' Check for environmental kills
            ' ... implementation
            
            Return KillMethod.Generic
        End Function
        
        Private Sub DropLoot()
            ' Drop weapons, trinkets, clues, etc.
            ' ... implementation
        End Sub
    End Class
    
    ' ============================================
    ' TIME LOOP MECHANICS
    ' ============================================
    
    Public Class LoopEventManager
        Private scheduledEvents As New Dictionary(Of Double, List(Of LoopEvent))()
        Private triggeredEvents As New List(Of LoopEvent)()
        
        Public Sub ScheduleEvent(eventTime As Double, loopEvent As LoopEvent)
            If Not scheduledEvents.ContainsKey(eventTime) Then
                scheduledEvents(eventTime) = New List(Of LoopEvent)()
            End If
            scheduledEvents(eventTime).Add(loopEvent)
        End Sub
        
        Public Sub Update(currentLoopTime As Double)
            ' Trigger events scheduled for this time
            For Each time In scheduledEvents.Keys
                If Math.Abs(time - currentLoopTime) < 0.1 Then ' Within 100ms
                    For Each loopEvent In scheduledEvents(time)
                        If Not loopEvent.HasTriggered OrElse loopEvent.RepeatsEachLoop Then
                            loopEvent.Trigger()
                            triggeredEvents.Add(loopEvent)
                        End If
                    Next
                End If
            Next
        End Sub
        
        Public Sub ResetForNewLoop()
            ' Reset triggered events for repeating events
            For Each loopEvent In triggeredEvents
                If loopEvent.RepeatsEachLoop Then
                    loopEvent.HasTriggered = False
                End If
            Next
            triggeredEvents.Clear()
        End Sub
    End Class
    
    Public Class LoopMemory
        Public Property RememberedLocations As New Dictionary(Of String, LocationMemory)()
        Public Property DiscoveredClues As New List(Of Clue)()
        Public Property LearnedPatterns As New List(Of Pattern)()
        Public Property LoopHistory As New List(Of LoopRecord)()
        
        Public Sub RecordLoop(loopRecord As LoopRecord)
            LoopHistory.Add(loopRecord)
            
            ' Keep only last 10 loops
            If LoopHistory.Count > 10 Then
                LoopHistory.RemoveAt(0)
            End If
        End Sub
        
        Public Sub RememberLocation(locationName As String, memory As LocationMemory)
            If RememberedLocations.ContainsKey(locationName) Then
                ' Update existing memory
                RememberedLocations(locationName).Merge(memory)
            Else
                RememberedLocations.Add(locationName, memory)
            End If
        End Sub
        
        Public Function GetLocationMemory(locationName As String) As LocationMemory
            If RememberedLocations.ContainsKey(locationName) Then
                Return RememberedLocations(locationName)
            End If
            Return New LocationMemory()
        End Function
    End Class
    
    ' ============================================
    ' ADVANCED PHYSICS EFFECTS
    ' ============================================
    
    Public Class TimeAnomaly
        Inherits PhysicsObject
        
        Public Property Type As TimeAnomalyType
        Public Property Strength As Single = 1.0F
        Public Property Radius As Single = 100.0F
        Public Property Duration As Double = 10.0 ' seconds
        
        Private activeTime As Double = 0
        
        Public Sub Update(deltaTime As Double)
            activeTime += deltaTime
            
            If activeTime >= Duration Then
                ' Anomaly expires
                RemoveFromWorld()
                Return
            End If
            
            ' Apply time distortion to nearby objects
            Dim nearbyObjects = GetObjectsInRadius(Radius)
            For Each obj In nearbyObjects
                ApplyTimeEffect(obj)
            Next
        End Sub
        
        Private Sub ApplyTimeEffect(obj As PhysicsObject)
            Select Case Type
                Case TimeAnomalyType.SlowField
                    obj.Velocity *= 0.5F
                    obj.AngularVelocity *= 0.5F
                    
                Case TimeAnomalyType.QuickField
                    obj.Velocity *= 1.5F
                    obj.AngularVelocity *= 1.5F
                    
                Case TimeAnomalyType.RewindField
                    ' Store and rewind object state
                    If TypeOf obj Is ITimeRewindable Then
                        CType(obj, ITimeRewindable).Rewind(Strength * 0.1F)
                    End If
                    
                Case TimeAnomalyType.StasisField
                    obj.Velocity = Vector2.Zero
                    obj.AngularVelocity = 0
            End Select
        End Sub
    End Class
    
    Public Class GravityWell
        Inherits PhysicsObject
        
        Public Property GravityStrength As Single = 500.0F
        Public Property Radius As Single = 200.0F
        Public Property IsBlackHole As Boolean = False
        Public Property PullStrength As Single = 1.0F
        
        Public Sub Update(deltaTime As Double)
            Dim nearbyObjects = GetObjectsInRadius(Radius)
            
            For Each obj In nearbyObjects
                If obj.IsAffectedByGravity Then
                    ApplyGravityPull(obj, deltaTime)
                End If
            Next
        End Sub
        
        Private Sub ApplyGravityPull(obj As PhysicsObject, deltaTime As Double)
            Dim direction = Vector2.Normalize(Position - obj.Position)
            Dim distance = Vector2.Distance(Position, obj.Position)
            
            If distance < 10 Then
                ' Object reached center - apply black hole effect
                If IsBlackHole Then
                    obj.Destroy()
                End If
                Return
            End If
            
            ' Inverse square law gravity
            Dim forceMagnitude = GravityStrength / (distance * distance)
            forceMagnitude *= PullStrength
            
            ' Apply force
            obj.Velocity += direction * forceMagnitude * deltaTime
            
            ' Apply tidal forces (stretch effect)
            If distance < Radius * 0.5 Then
                ApplyTidalForces(obj, direction, distance)
            End If
        End Sub
        
        Private Sub ApplyTidalForces(obj As PhysicsObject, direction As Vector2, distance As Single)
            ' Stretch object along radial direction
            Dim stretchFactor = 1.0F + (Radius * 0.5F - distance) / 100.0F
            obj.Scale = New Vector2(stretchFactor, 1.0F / stretchFactor)
            obj.Rotation = Math.Atan2(direction.Y, direction.X)
        End Sub
    End Class
    
    Public Class DestructibleObject
        Inherits PhysicsObject
        
        Public Property Health As Single = 100.0F
        Public Property DestructionThreshold As Single = 25.0F
        Public Property IsDestroyed As Boolean = False
        Public Property DebrisCount As Integer = 10
        Public Property ExplosionForce As Single = 1000.0F
        Public Property CanRegenerate As Boolean = False
        Public Property RegenerationTime As Double = 30.0
        
        Private regenerationTimer As Double = 0
        
        Public Overrides Sub OnCollision(other As PhysicsObject, collisionInfo As CollisionInfo)
            MyBase.OnCollision(other, collisionInfo)
            
            ' Calculate collision damage
            Dim relativeVelocity = other.Velocity - Me.Velocity
            Dim impactForce = Vector2.Dot(relativeVelocity, collisionInfo.Normal) * Me.Mass
            
            If impactForce > DestructionThreshold Then
                TakeDamage(impactForce)
            End If
        End Sub
        
        Public Sub TakeDamage(damage As Single)
            Health -= damage
            
            If Health <= 0 AndAlso Not IsDestroyed Then
                Destroy()
            End If
        End Sub
        
        Public Sub Destroy()
            IsDestroyed = True
            
            ' Create explosion effect
            CreateExplosion()
            
            ' Spawn debris
            CreateDebris()
            
            ' Trigger world event
            worldEvents.Trigger("ObjectDestroyed", Me)
        End Sub
        
        Private Sub CreateExplosion()
            ' Apply radial force to nearby objects
            Dim nearbyObjects = GetObjectsInRadius(200)
            
            For Each obj In nearbyObjects
                Dim direction = Vector2.Normalize(obj.Position - Position)
                Dim distance = Vector2.Distance(obj.Position, Position)
                Dim force = ExplosionForce / (distance + 1)
                
                obj.Velocity += direction * force
                obj.TakeDamage(force * 0.1F, DamageType.Explosive, Nothing)
            Next
            
            ' Create visual explosion
            ' ... particle system
        End Sub
        
        Private Sub CreateDebris()
            For i As Integer = 1 To DebrisCount
                Dim debris As New PhysicsObject() With {
                    .Position = Position + New Vector2(Rnd() * 50 - 25, Rnd() * 50 - 25),
                    .Velocity = New Vector2(Rnd() * 200 - 100, Rnd() * 200 - 100),
                    .Size = New Vector2(5, 5),
                    .Mass = 1.0F,
                    .CollisionLayer = CollisionLayer.Debris
                }
                
                ' Add debris to physics world
                AddPhysicsObject(debris)
            Next
        End Sub
        
        Public Sub Update(deltaTime As Double)
            If IsDestroyed AndAlso CanRegenerate Then
                regenerationTimer += deltaTime
                
                If regenerationTimer >= RegenerationTime Then
                    Regenerate()
                End If
            End If
            
            MyBase.Update(deltaTime)
        End Sub
        
        Private Sub Regenerate()
            IsDestroyed = False
            Health = 100.0F
            regenerationTimer = 0
            Opacity = 1.0F
            worldEvents.Trigger("ObjectRegenerated", Me)
        End Sub
    End Class
    
    Public Class RagdollPhysics
        Inherits PhysicsObject
        
        Public Property Bones As New List(Of RagdollBone)()
        Public Property Joints As New List(Of PhysicsJoint)()
        Public Property IsActive As Boolean = True
        Public Property DecayTimer As Double = 30.0 ' Seconds until removal
        
        Public Sub New(baseObject As GameObject)
            ' Create bones from object's skeleton
            CreateBonesFromSkeleton(baseObject)
            
            ' Create joints between bones
            CreateJoints()
            
            ' Apply initial forces
            ApplyInitialForces(baseObject)
        End Sub
        
        Private Sub CreateBonesFromSkeleton(baseObject As GameObject)
            ' Simplified: Create main body parts
            Bones.Add(New RagdollBone("Head", baseObject.Position + New Vector2(0, -20)))
            Bones.Add(New RagdollBone("Torso", baseObject.Position))
            Bones.Add(New RagdollBone("LeftArm", baseObject.Position + New Vector2(-30, 0)))
            Bones.Add(New RagdollBone("RightArm", baseObject.Position + New Vector2(30, 0)))
            Bones.Add(New RagdollBone("LeftLeg", baseObject.Position + New Vector2(-10, 30)))
            Bones.Add(New RagdollBone("RightLeg", baseObject.Position + New Vector2(10, 30)))
        End Sub
        
        Private Sub CreateJoints()
            ' Connect bones with constraints
            Joints.Add(New DistanceJoint(Bones(0), Bones(1), 20)) ' Head to torso
            Joints.Add(New DistanceJoint(Bones(1), Bones(2), 25)) ' Torso to left arm
            Joints.Add(New DistanceJoint(Bones(1), Bones(3), 25)) ' Torso to right arm
            Joints.Add(New DistanceJoint(Bones(1), Bones(4), 30)) ' Torso to left leg
            Joints.Add(New DistanceJoint(Bones(1), Bones(5), 30)) ' Torso to right leg
        End Sub
        
        Private Sub ApplyInitialForces(baseObject As GameObject)
            ' Transfer momentum from original object
            For Each bone In Bones
                bone.Velocity = baseObject.Velocity + 
                    New Vector2(Rnd() * 100 - 50, Rnd() * 100 - 50)
                bone.AngularVelocity = Rnd() * 10 - 5
            Next
        End Sub
        
        Public Overrides Sub Update(deltaTime As Double)
            If Not IsActive Then Return
            
            ' Update bones
            For Each bone In Bones
                bone.Update(deltaTime)
            Next
            
            ' Solve joints
            For i As Integer = 1 To 3 ' Multiple solver iterations
                For Each joint In Joints
                    joint.Solve()
                Next
            Next
            
            ' Apply decay
            DecayTimer -= deltaTime
            If DecayTimer <= 0 Then
                IsActive = False
            End If
        End Sub
    End Class
    
    ' ============================================
    ' GAME LOOP AND RENDERING
    ' ============================================
    
    Private Sub DeathLoopGame_Load(sender As Object, e As EventArgs) Handles MyBase.Load
        ' Initialize game
        InitializeGame()
        
        ' Set up timers
        gameTimer.Interval = CInt(RENDER_TIME_STEP * 1000)
        physicsTimer.Interval = CInt(PHYSICS_TIME_STEP * 1000)
        
        gameTimer.Start()
        physicsTimer.Start()
        
        ' Set up rendering
        renderTarget = New Bitmap(ClientSize.Width, ClientSize.Height)
        renderGraphics = Graphics.FromImage(renderTarget)
        renderGraphics.SmoothingMode = SmoothingMode.AntiAlias
        renderGraphics.InterpolationMode = InterpolationMode.HighQualityBicubic
        
        ' Enable double buffering
        Me.DoubleBuffered = True
    End Sub
    
    Private Sub InitializeGame()
        ' Initialize player
        player = New DeathLoopPlayer() With {
            .Position = New Vector2(100, 100),
            .Size = New Vector2(30, 60),
            .Mass = 70.0F ' 70kg
        }
        
        ' Initialize visionaries
        InitializeVisionaries()
        
        ' Initialize world events
        InitializeWorldEvents()
        
        ' Initialize physics world
        physicsWorld = New PhysicsWorld()
        
        ' Set up initial loop
        ResetLoop()
    End Sub
    
    Private Sub InitializeVisionaries()
        ' Create 8 visionaries (as in DeathLoop)
        For i As Integer = 1 To 8
            Dim visionary As New Visionary() With {
                .Name = $"Visionary_{i}",
                .Position = New Vector2(200 + i * 100, 200),
                .Size = New Vector2(40, 80),
                .Mass = 80.0F
            }
            
            ' Set up unique abilities based on visionary type
            visionary.UniqueAbility = CreateVisionaryAbility(i)
            
            visionaries.Add(visionary)
        Next
    End Sub
    
    Private Function CreateVisionaryAbility(visionaryIndex As Integer) As VisionaryAbility
        Select Case visionaryIndex
            Case 1 ' Julianna - Time Manipulation
                Return New TimeManipulationAbility()
            Case 2 ' Aleksis - Strength Enhancement
                Return New StrengthEnhancementAbility()
            Case 3 ' Harriet - Poison Cloud
                Return New PoisonCloudAbility()
            ' ... other visionaries
            Case Else
                Return New GenericAbility()
        End Select
    End Function
    
    Private Sub InitializeWorldEvents()
        ' Schedule daily events
        loopEvents.ScheduleEvent(3600, New LoopEvent("MorningAnnouncement", 
            Sub() worldEvents.Trigger("MorningAnnouncement", "Good morning, Blackreef!")))
        
        loopEvents.ScheduleEvent(14400, New LoopEvent("LunchTime", 
            Sub()
                For Each v In visionaries
                    If v.IsAlive Then
                        v.Routine.SetCurrentActivity("Eating")
                    End If
                Next
            End Sub))
        
        loopEvents.ScheduleEvent(64800, New LoopEvent("EveningParty", 
            Sub() worldEvents.Trigger("EveningParty", "The party at Aleksis' mansion begins!")))
    End Sub
    
    Private Sub gameTimer_Tick(sender As Object, e As EventArgs) Handles gameTimer.Tick
        ' Update game logic
        UpdateGame(RENDER_TIME_STEP)
        
        ' Render
        RenderGame()
        
        ' Invalidate to trigger paint
        Invalidate()
    End Sub
    
    Private Sub physicsTimer_Tick(sender As Object, e As EventArgs) Handles physicsTimer.Tick
        ' Update physics
        UpdatePhysics(PHYSICS_TIME_STEP)
    End Sub
    
    Private Sub UpdateGame(deltaTime As Double)
        ' Update game time
        gameTime += deltaTime
        loopTime += deltaTime
        
        ' Check for loop reset
        If loopTime >= LOOP_DURATION Then
            ResetLoop()
            Return
        End If
        
        ' Update player
        Dim input = GetPlayerInput()
        player.Update(deltaTime, input)
        
        ' Update visionaries
        For Each visionary In visionaries
            If visionary.IsAlive Then
                visionary.Update(deltaTime, player.Position)
            End If
        Next
        
        ' Update time anomalies
        For Each anomaly In timeAnomalies
            anomaly.Update(deltaTime)
        Next
        
        ' Update gravity wells
        For Each gravityWell In gravityWells
            gravityWell.Update(deltaTime)
        Next
        
        ' Update destructible objects
        For Each obj In destructibleObjects
            obj.Update(deltaTime)
        Next
        
        ' Update ragdolls
        For Each ragdoll In ragdolls
            ragdoll.Update(deltaTime)
        Next
        
        ' Update loop events
        loopEvents.Update(loopTime)
        
        ' Update world events
        worldEvents.Update(deltaTime)
        
        ' Check win condition
        CheckWinCondition()
    End Sub
    
    Private Sub UpdatePhysics(deltaTime As Double)
        physicsWorld.Update(deltaTime)
    End Sub
    
    Private Sub RenderGame()
        ' Clear background
        renderGraphics.Clear(Color.FromArgb(20, 20, 30))
        
        ' Draw game world
        DrawWorld(renderGraphics)
        
        ' Draw physics objects
        DrawPhysicsObjects(renderGraphics)
        
        ' Draw player
        DrawPlayer(renderGraphics)
        
        ' Draw visionaries
        DrawVisionaries(renderGraphics)
        
        ' Draw UI
        DrawUI(renderGraphics)
        
        ' Draw time loop overlay
        DrawTimeLoopOverlay(renderGraphics)
    End Sub
    
    Private Sub DrawWorld(g As Graphics)
        ' Draw terrain
        Using terrainBrush As New SolidBrush(Color.FromArgb(40, 40, 60))
            g.FillRectangle(terrainBrush, 0, 0, ClientSize.Width, ClientSize.Height)
        End Using
        
        ' Draw buildings and structures
        ' ... implementation
        
        ' Draw time anomalies
        For Each anomaly In timeAnomalies
            DrawTimeAnomaly(g, anomaly)
        Next
        
        ' Draw gravity wells
        For Each gravityWell In gravityWells
            DrawGravityWell(g, gravityWell)
        Next
    End Sub
    
    Private Sub DrawPlayer(g As Graphics)
        ' Draw player body
        Using playerBrush As New SolidBrush(Color.FromArgb(100, 150, 255))
            g.FillRectangle(playerBrush,
                           player.Position.X - player.Size.X / 2,
                           player.Position.Y - player.Size.Y / 2,
                           player.Size.X, player.Size.Y)
        End Using
        
        ' Draw player health bar
        DrawHealthBar(g, player.Position, player.Health, player.MaxHealth, Color.Green)
        
        ' Draw player energy bar
        DrawEnergyBar(g, player.Position, player.Energy, player.MaxEnergy, Color.Blue)
        
        ' Draw remembered path
        If player.RememberedPaths.Count > 1 Then
            Using pathPen As New Pen(Color.FromArgb(100, 255, 255, 255), 2)
                For i As Integer = 0 To player.RememberedPaths.Count - 2
                    g.DrawLine(pathPen,
                              player.RememberedPaths(i),
                              player.RememberedPaths(i + 1))
                Next
            End Using
        End If
        
        ' Draw echo if active
        If player.CanCreateEcho AndAlso player.IsRecordingEcho Then
            Using echoPen As New Pen(Color.FromArgb(150, 255, 255, 0), 1)
                For i As Integer = 0 To player.echoRecording.Count - 2
                    Dim frame1 = player.echoRecording(i)
                    Dim frame2 = player.echoRecording(i + 1)
                    g.DrawLine(echoPen, frame1.Position, frame2.Position)
                Next
            End Using
        End If
    End Sub
    
    Private Sub DrawVisionaries(g As Graphics)
        For Each visionary In visionaries
            If visionary.IsAlive Then
                ' Draw visionary
                Using visionaryBrush As New SolidBrush(Color.FromArgb(255, 100, 100))
                    g.FillRectangle(visionaryBrush,
                                   visionary.Position.X - visionary.Size.X / 2,
                                   visionary.Position.Y - visionary.Size.Y / 2,
                                   visionary.Size.X, visionary.Size.Y)
                End Using
                
                ' Draw name
                Using font As New Font("Arial", 8)
                    g.DrawString(visionary.Name, font, Brushes.White,
                                visionary.Position.X - 30,
                                visionary.Position.Y - 40)
                End Using
                
                ' Draw awareness bar
                DrawAwarenessBar(g, visionary.Position, visionary.Awareness)
                
                ' Draw health bar
                DrawHealthBar(g, visionary.Position, visionary.Health, 
                             visionary.MaxHealth, Color.Red)
            End If
        Next
    End Sub
    
    Private Sub DrawUI(g As Graphics)
        ' Draw loop timer
        Dim timeLeft = LOOP_DURATION - loopTime
        Dim timeString = $"Loop Time: {FormatTime(timeLeft)}"
        
        Using font As New Font("Arial", 12, FontStyle.Bold)
            g.DrawString(timeString, font, Brushes.White, 10, 10)
            g.DrawString($"Loop: {loopCount}", font, Brushes.Yellow, 10, 30)
            g.DrawString($"Deaths: {player.LoopDeaths}", font, Brushes.Red, 10, 50)
            g.DrawString($"Kills: {player.LoopKills}", font, Brushes.Green, 10, 70)
        End Using
        
        ' Draw objectives
        Dim aliveCount = visionaries.Count(Function(v) v.IsAlive)
        g.DrawString($"Visionaries Remaining: {aliveCount}", 
                    New Font("Arial", 10), Brushes.Orange, 10, 100)
        
        ' Draw discovered clues
        If player.DiscoveredSecrets.Count > 0 Then
            g.DrawString($"Secrets Found: {player.DiscoveredSecrets.Count}", 
                        New Font("Arial", 10), Brushes.Cyan, 10, 120)
        End If
    End Sub
    
    Private Sub DrawTimeLoopOverlay(g As Graphics)
        ' Draw temporal distortion effects
        If physicsWorld.TimeDilation <> 1.0F Then
            Using overlayBrush As New SolidBrush(
                Color.FromArgb(50, 
                    If(physicsWorld.TimeDilation < 1, 0, 255),
                    If(physicsWorld.TimeDilation < 1, 0, 150),
                    255))
                g.FillRectangle(overlayBrush, 0, 0, ClientSize.Width, ClientSize.Height)
            End Using
            
            ' Draw time dilation indicator
            Using font As New Font("Arial", 14, FontStyle.Bold)
                Dim timeText = If(physicsWorld.TimeDilation < 0, 
                                 $"TIME REWIND: {Math.Abs(physicsWorld.TimeDilation):P0}",
                                 $"TIME DILATION: {physicsWorld.TimeDilation:P0}")
                
                g.DrawString(timeText, font, Brushes.White, 
                            ClientSize.Width / 2 - 100, 20)
            End Using
        End If
    End Sub
    
    Private Function GetPlayerInput() As PlayerInput
        ' Get keyboard state
        Dim state = Keyboard.GetState()
        
        Return New PlayerInput() With {
            .MoveUp = state.IsKeyDown(Keys.W) Or state.IsKeyDown(Keys.Up),
            .MoveDown = state.IsKeyDown(Keys.S) Or state.IsKeyDown(Keys.Down),
            .MoveLeft = state.IsKeyDown(Keys.A) Or state.IsKeyDown(Keys.Left),
            .MoveRight = state.IsKeyDown(Keys.D) Or state.IsKeyDown(Keys.Right),
            .IsSprinting = state.IsKeyDown(Keys.LeftShift),
            .IsCrouching = state.IsKeyDown(Keys.LeftControl),
            .IsJumping = state.IsKeyDown(Keys.Space),
            .IsShooting = state.IsKeyDown(Keys.MouseLeft),
            .IsMeleeAttacking = state.IsKeyDown(Keys.F),
            .IsRewindingTime = state.IsKeyDown(Keys.Q),
            .IsSlowingTime = state.IsKeyDown(Keys.E),
            .IsFreezingTime = state.IsKeyDown(Keys.R),
            .IsPhasing = state.IsKeyDown(Keys.C),
            .IsRecordingEcho = state.IsKeyDown(Keys.V),
            .SwitchWeaponForward = state.IsKeyDown(Keys.Tab),
            .SwitchWeaponBackward = state.IsKeyDown(Keys.LeftAlt)
        }
    End Function
    
    Private Sub ResetLoop()
        loopCount += 1
        loopTime = 0
        
        ' Reset visionaries (except those permanently dead in story mode)
        For Each visionary In visionaries
            If visionary.IsAlive = False Then
                ' Check if visionary should respawn
                If loopCount Mod 3 <> 0 Then ' Every 3 loops
                    visionary.IsAlive = True
                    visionary.Health = visionary.MaxHealth
                    visionary.Awareness = 0
                    visionary.CurrentState = AIState.Patrolling
                End If
            End If
        Next
        
        ' Reset player position but keep memories
        player.Position = New Vector2(100, 100)
        player.Health = player.MaxHealth
        player.Energy = player.MaxEnergy
        player.IsInvulnerable = False
        
        ' Reset physics objects
        destructibleObjects.RemoveAll(Function(o) o.IsDestroyed AndAlso Not o.CanRegenerate)
        
        ' Regenerate objects that can
        For Each obj In destructibleObjects
            If obj.CanRegenerate Then
                obj.Regenerate()
            End If
        Next
        
        ' Reset loop events
        loopEvents.ResetForNewLoop()
        
        ' Record loop in memory
        playerMemory.RecordLoop(New LoopRecord(loopCount, player.LoopKills, 
                                              player.LoopDeaths, loopTime))
        
        gameState = GameState.LoopReset
    End Sub
    
    Private Sub CheckWinCondition()
        ' Check if all visionaries are dead
        If visionaries.All(Function(v) Not v.IsAlive) Then
            gameState = GameState.VisionaryKilled
            ' Trigger ending sequence
            TriggerEndingSequence()
        End If
    End Sub
    
    Private Sub TriggerEndingSequence()
        ' Freeze time
        physicsWorld.TimeDilation = 0.001F
        
        ' Show ending cinematic
        ' ... implementation
    End Sub
    
    Private Function FormatTime(seconds As Double) As String
        Dim hours = CInt(Math.Floor(seconds / 3600))
        Dim minutes = CInt(Math.Floor((seconds Mod 3600) / 60))
        Dim secs = CInt(Math.Floor(seconds Mod 60))
        
        Return $"{hours:D2}:{minutes:D2}:{secs:D2}"
    End Function
    
    ' Helper drawing methods
    Private Sub DrawHealthBar(g As Graphics, position As Vector2, 
                             currentHealth As Single, maxHealth As Single, color As Color)
        Dim barWidth As Integer = 40
        Dim barHeight As Integer = 5
        Dim barX As Single = position.X - barWidth / 2
        Dim barY As Single = position.Y - 50
        
        ' Draw background
        g.FillRectangle(Brushes.DarkRed, barX, barY, barWidth, barHeight)
        
        ' Draw health
        Dim healthPercent = currentHealth / maxHealth
        g.FillRectangle(New SolidBrush(color), 
                       barX, barY, barWidth * healthPercent, barHeight)
        
        ' Draw border
        g.DrawRectangle(Pens.Black, barX, barY, barWidth, barHeight)
    End Sub
    
    Private Sub DrawEnergyBar(g As Graphics, position As Vector2, 
                             currentEnergy As Single, maxEnergy As Single, color As Color)
        Dim barWidth As Integer = 40
        Dim barHeight As Integer = 3
        Dim barX As Single = position.X - barWidth / 2
        Dim barY As Single = position.Y - 45
        
        ' Draw background
        g.FillRectangle(Brushes.DarkBlue, barX, barY, barWidth, barHeight)
        
        ' Draw energy
        Dim energyPercent = currentEnergy / maxEnergy
        g.FillRectangle(New SolidBrush(color), 
                       barX, barY, barWidth * energyPercent, barHeight)
    End Sub
    
    Private Sub DrawAwarenessBar(g As Graphics, position As Vector2, awareness As Single)
        Dim barWidth As Integer = 40
        Dim barHeight As Integer = 3
        Dim barX As Single = position.X - barWidth / 2
        Dim barY As Single = position.Y - 60
        
        ' Draw awareness
        Dim awarenessColor = Color.FromArgb(
            CInt(awareness * 2.55),
            CInt(awareness * 2.55),
            0)
        
        g.FillRectangle(New SolidBrush(awarenessColor), 
                       barX, barY, barWidth * awareness / 100, barHeight)
    End Sub
    
    Protected Overrides Sub OnPaint(e As PaintEventArgs)
        MyBase.OnPaint(e)
        
        ' Draw the rendered frame
        If renderTarget IsNot Nothing Then
            e.Graphics.DrawImage(renderTarget, 0, 0)
        End If
    End Sub
    
    Protected Overrides Sub OnKeyDown(e As KeyEventArgs)
        MyBase.OnKeyDown(e)
        
        ' Handle debug keys
        Select Case e.KeyCode
            Case Keys.F1
                ' Toggle physics debug
            Case Keys.F2
                ' Spawn time anomaly
                SpawnTimeAnomaly(TimeAnomalyType.SlowField, player.Position)
            Case Keys.F3
                ' Spawn gravity well
                SpawnGravityWell(player.Position)
            Case Keys.F4
                ' Kill all visionaries (debug)
                For Each v In visionaries
                    v.Health = 0
                Next
            Case Keys.F5
                ' Reset loop manually
                ResetLoop()
            Case Keys.P
                ' Pause game
                gameState = If(gameState = GameState.Running, 
                              GameState.Paused, GameState.Running)
        End Select
    End Sub
    
    Private Sub SpawnTimeAnomaly(type As TimeAnomalyType, position As Vector2)
        Dim anomaly As New TimeAnomaly() With {
            .Type = type,
            .Position = position,
            .Strength = 1.0F,
            .Radius = 150.0F,
            .Duration = 15.0
        }
        timeAnomalies.Add(anomaly)
    End Sub
    
    Private Sub SpawnGravityWell(position As Vector2)
        Dim gravityWell As New GravityWell() With {
            .Position = position,
            .GravityStrength = 800.0F,
            .Radius = 250.0F,
            .IsBlackHole = False,
            .PullStrength = 1.0F
        }
        gravityWells.Add(gravityWell)
    End Sub
End Class

' ============================================
' SUPPORTING CLASSES
' ============================================

Public Class PlayerInput
    Public Property MoveUp As Boolean
    Public Property MoveDown As Boolean
    Public Property MoveLeft As Boolean
    Public Property MoveRight As Boolean
    Public Property IsSprinting As Boolean
    Public Property IsCrouching As Boolean
    Public Property IsJumping As Boolean
    Public Property IsShooting As Boolean
    Public Property IsMeleeAttacking As Boolean
    Public Property IsRewindingTime As Boolean
    Public Property IsSlowingTime As Boolean
    Public Property IsFreezingTime As Boolean
    Public Property IsPhasing As Boolean
    Public Property IsRecordingEcho As Boolean
    Public Property SwitchWeaponForward As Boolean
    Public Property SwitchWeaponBackward As Boolean
    Public Property AimDirection As Vector2
End Class

Public Enum TimeAnomalyType
    SlowField
    QuickField
    RewindField
    StasisField
End Enum

Public Enum DamageType
    Kinetic
    Energy
    Explosive
    Poison
    Fire
    Ice
    Electric
    Void
End Enum

Public Enum AIState
    Patrolling
    Investigating
    Combat
    Fleeing
    UsingAbility
    Dead
End Enum

Public Enum KillMethod
    None
    Generic
    WeaponSpecific
    Environmental
    Stealth
    TimeLoop
End Enum

Public Class LoopEvent
    Public Property Name As String
    Public Property TriggerTime As Double
    Public Property Action As Action
    Public Property HasTriggered As Boolean = False
    Public Property RepeatsEachLoop As Boolean = False
    
    Public Sub New(name As String, action As Action)
        Me.Name = name
        Me.Action = action
    End Sub
    
    Public Sub Trigger()
        Action.Invoke()
        HasTriggered = True
    End Sub
End Class

Public Class LocationMemory
    Public Property DiscoveredItems As New List(Of String)()
    Public Property KnownPaths As New List(Of Vector2)()
    Public Property HiddenAreas As New List(Of String)()
    Public Property DangerZones As New List(Of String)()
    
    Public Sub Merge(other As LocationMemory)
        DiscoveredItems.AddRange(other.DiscoveredItems)
        KnownPaths.AddRange(other.KnownPaths)
        HiddenAreas.AddRange(other.HiddenAreas)
        DangerZones.AddRange(other.DangerZones)
    End Sub
End Class

Public Class LoopRecord
    Public Property LoopNumber As Integer
    Public Property Kills As Integer
    Public Property Deaths As Integer
    Public Property TimeSurvived As Double
    Public Property DateRecorded As DateTime
    
    Public Sub New(loopNumber As Integer, kills As Integer, 
                  deaths As Integer, timeSurvived As Double)
        Me.LoopNumber = loopNumber
        Me.Kills = kills
        Me.Deaths = deaths
        Me.TimeSurvived = timeSurvived
        Me.DateRecorded = DateTime.Now
    End Sub
End Class

Public Class PlayerFrame
    Public Property Position As Vector2
    Public Property Velocity As Vector2
    Public Property Rotation As Single
    Public Property Action As String
    Public Property Timestamp As Double
End Class

' Entry point
Module Program
    <STAThread>
    Sub Main()
        Application.EnableVisualStyles()
        Application.SetCompatibleTextRenderingDefault(False)
        Application.Run(New DeathLoopGame())
    End Sub
End Module