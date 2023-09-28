import Icon from './svg/bone.svg';

function BoneIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="bone icon" />
    </div>
  );
}

export default BoneIcon;
